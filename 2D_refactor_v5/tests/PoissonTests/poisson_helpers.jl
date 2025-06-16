"""
solve on flat periodic domain using FE space with zero mean constraint
"""
function solve_poisson_periodic(domain,partition,p,degree,uex)

  model = CartesianDiscreteModel(domain,partition,isperiodic=ntuple(x->true,length(partition)))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  ## assess zero mean of analytic function
  println("Zero mean is: " , sum(∫(uex)dΩ))

  #### FE Problem

  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,
                        constraint=:zeromean)
  U = TrialFESpace(V)

  rhs(x) = -1.0*(laplacian(uex)(x))

  poisson_biform(u,v) = ∫( gradient(u)⋅gradient(v) )dΩ
  poisson_liform(v) =  ∫( rhs*v )dΩ

  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

  b = get_vector(op)
  println("Compatibility: ", sum(b))

  uh = solve(LUSolver(),op)

  e = sum(∫((uh-uex)⊙(uh-uex))dΩ)

  println("Errors: ", e)
  return e
end

"""
solve on flat periodic domain using larange multipliers to enforce zeromean constraint
"""
function solve_poisson_periodic_lagrange(domain,partition,p,degree,uex)

  model = CartesianDiscreteModel(domain,partition,isperiodic=ntuple(x->true,length(partition)))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  ## assess zero mean of analytic function
  println("Zero mean is: " , sum(∫(uex)dΩ))

  #### FE Problem

  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
  U = TrialFESpace(V)
  Λ = ConstantFESpace(model)
  M = TrialFESpace(Λ)
  X = MultiFieldFESpace([U,M])
  Y = MultiFieldFESpace([V,Λ])

  rhs(x) = -1.0*(laplacian(uex)(x))

  poisson_biform((u,μ),(v,λ)) = ∫( gradient(u)⋅gradient(v)  )dΩ  + ∫(v*μ)dΩ + ∫(λ*u)dΩ
  poisson_liform((v,λ)) = ∫( rhs*v )dΩ #+ ∫(λ*uex)dΩ

  op = AffineFEOperator(poisson_biform,poisson_liform,X,Y)

  b = get_vector(op)
  println("Compatibility: ", sum(b))

  uh,μh = solve(LUSolver(),op)

  e = sum(∫((uh-uex)⊙(uh-uex))dΩ)

  println("Errors: ", e)
  return e
end


function poisson_periodic(domain,partition,p,degree,uex;lagrange=false,uzeromean=false)
  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,partition,isperiodic=ntuple(x->true,length(partition))))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  rhs(x) = -1.0*(laplacian(uex)(x))

  ### FE problem - multiifield
  if lagrange
    V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
    U = TrialFESpace(V)
    Λ = ConstantFESpace(model)
    M = TrialFESpace(Λ)
    X = MultiFieldFESpace([U,M])
    Y = MultiFieldFESpace([V,Λ])
    poisson_biformX((u,μ),(v,λ)) = ∫( gradient(u)⋅gradient(v)  )dΩ  + ∫(v*μ)dΩ + ∫(λ*u)dΩ

    function poisson_liformY((v,λ))
      if uzeromean # force uex to have zeromean
        return ∫( rhs*v )dΩ  + ∫(λ*uex)dΩ
      else # only force u to have zero mean
        return  ∫( rhs*v )dΩ
      end
    end

    op = AffineFEOperator(poisson_biformX,poisson_liformY,X,Y)
    uh,μh = solve(BackslashSolver(),op)

    # b = get_matrix(op)
    # println("Compatibility: ", b)

    return sum(∫((uh-uex)⊙(uh-uex))dΩ)
  end


  ### FE problem -- single field
  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,constraint=:zeromean)
  U = TrialFESpace(V)
  poisson_biform(u,v) = ∫( gradient(u)⋅gradient(v) )dΩ
  poisson_liform(v) =  ∫( rhs*v )dΩ
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
  uh = solve(BackslashSolver(),op)
  return sum(∫((uh-uex)⊙(uh-uex))dΩ)
end



"""
solve on manifold using surface diff operators
"""
function solve_poisson_manifold(domain,n,p,degree,u,metric_func;isperiodic=ntuple(i->false,length(n)))

  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,n,isperiodic=isperiodic))

  Ω = Triangulation(model)
  m = Metric(metric_func,Ω)

  dΩ = Measure(Ω,degree)
  dΩg =  Measure(m,Ω,degree)

  #### FE Problem
  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,
                        dirichlet_tags="boundary")
  U = TrialFESpace(V,u)

  if isperiodic == (true,true)
    println("zero mean")
    V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,constraint=:zeromean)
    U = TrialFESpace(V)
  end

  ucf = CellField(u,Ω)

  rhs = -1.0*surface_laplacian(ucf,m)

  poisson_biform(u,v) = ∫( surface_gradient(u,m)⋅gradient(v) )dΩg
  poisson_liform(v) =  ∫( (rhs*v) )dΩg

  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

  b = get_vector(op)
  println("Compatibility: ", sum(b))

  uh = solve(LUSolver(),op)

  e = l2(uh-ucf,dΩ)
  eg = l2(uh-ucf,dΩg)

  println("Errors: ", e, "; ", eg)
  return e, eg
end




function plot_error(ns,errs,leginf;
  ls=[:solid, :dash, :dot, :dashdot, :dashdotdot],
  colors = palette(:tab10),
  markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross] )

  nsims = Int(length(errs)/length(ns))
  for i in 1:nsims
    idx1 = 1 + (i-1)*length(ns)
    idx2 = (i)*length(ns)
    plot!(ns,
          errs[idx1:idx2],
          lw=3,
          c=colors[i],ls=ls[i], markershape=markers[i],
          label=leginf[i])
  end

end



function plot_convergence(ns,errs,errs_g,leginf)
  plot()
  plot_error(ns,errs,leginf;
    ls = fill(:solid,length(leginf)),
    markers = fill(:circle,length(leginf)))

  plot_error(ns,errs_g,fill(false,length(leginf));
    ls = fill(:dash,length(leginf)),
    markers = fill(:xcross,length(leginf)))

  plot!(yscale=:log10,framestyle=:box,
  xscale=:log10,
  xlabel="n cells",
  ylabel=L"L2(u_{ex} - u_h)"
  )
  plot!(show=true)
  plot!(xtickfontsize=10,ytickfontsize=10,
  legendfontsize=10,guidefontsize=10)
end
