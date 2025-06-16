"""
solve on poisson equation on flat periodic domain, doubly periodic
  * assume u_ex satisfies the compatibility condition ∫ Δuex = 0


By default, solve for u ∈ FESpace(...;constraint=:zeromean)
  * this method assumes u_ex is already zero mean

kwargs
  - lagrange == true  -> use lagrange multiplier to enforce ∫u = 0
      * this method assumes u_ex is already zero mean
  - uzeromean == true -> use lagrange multipler to enforce ∫uex = 0
"""
function solve_poisson_periodic(domain,partition,p,degree,uex;lagrange=false,uzeromean=false)
  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,partition,isperiodic=ntuple(x->true,length(partition))))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  rhs(x) = -1.0*(laplacian(uex)(x))

  ### FE problem - multiifield with lagrange multipliers
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
    uh,μh = solve(LUSolver(),op)

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
  uh = solve(LUSolver(),op)
  return sum(∫((uh-uex)⊙(uh-uex))dΩ)
end

"""
solve on poisson equation on flat periodic domain, doubly periodic, in mixed form
  σ = ∇ ⋅ u
  ∇⋅σ = -f
where f = -Δuex. Assume uex satisfies periodic boundary conditions

This method uses lagrange mulitplers to
  1. enforce the compatibility condition ∫ Δuex = 0
  2. enforce the zero means condition ∫ uex = ∫ u = 0
"""
function solve_poisson_dual_form(domain,partition,p,degree,uex)

  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,partition,isperiodic=ntuple(x->true,length(partition))))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)


  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:L2)
  U = TrialFESpace(V)

  T = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p); conformity=:Hdiv)
  S = TrialFESpace(T)

  Λ = ConstantFESpace(model)
  M = TrialFESpace(Λ)

  X = MultiFieldFESpace([S,U,M])
  Y = MultiFieldFESpace([T,V,Λ])

  ### force compatibility
  sigma_ex(x) = gradient(uex)(x)

  _X = MultiFieldFESpace([S,M])
  _Y = MultiFieldFESpace([T,Λ])

  biformS((s,μ),(t,λ)) = ∫( s⋅t )dΩ + ∫( divergence(s)*λ )dΩ  + ∫( divergence(t)*μ )dΩ
  liformS((t,λ)) = ∫( sigma_ex ⋅ t )dΩ
  op = AffineFEOperator(biformS,liformS,_X,_Y)
  sigma_exh,μh = solve(LUSolver(),op)

  # sum(∫((sigma_exh-sigma_ex)⊙(sigma_exh-sigma_ex))dΩ)


  ### # dual form -- with periodicity, force zero mean
  rhs = -1.0*divergence(sigma_exh)

  biformX((s,u,μ),(t,v,λ)) = (  ∫( s⋅t + divergence(t)*u )dΩ
                            + ∫( divergence(s)*v  )dΩ
                            + ∫(v*μ)dΩ + ∫(λ*u)dΩ
                      )
  liformY((t,v,λ)) = ∫( -(rhs*v) )dΩ  + ∫(λ*uex)dΩ

  op = AffineFEOperator(biformX,liformY,X,Y)
  sh,uh,μh = solve(LUSolver(),op)

  e = sum(∫((uh-uex)⊙(uh-uex))dΩ)
  println("Error = ", e)
  return e

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




function plot_error(ns,errs;
  leginf = fill(false,Int(length(errs)/length(ns))),
  ls=[:solid, :dash, :dot, :dashdot, :dashdotdot],
  colors = palette(:tab10),
  markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross],
  ms=[6,6,8,6,6,8,8] )

  nsims = Int(length(errs)/length(ns))
  for i in 1:nsims
    idx1 = 1 + (i-1)*length(ns)
    idx2 = (i)*length(ns)
    plot!(ns,
          errs[idx1:idx2],
          lw=3,
          markersize=ms[i],
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
