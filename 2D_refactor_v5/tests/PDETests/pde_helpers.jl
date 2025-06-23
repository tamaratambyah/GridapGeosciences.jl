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
function solve_poisson_periodic_dual_form(domain,partition,p,degree,uex)

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

Note this method is NOT applicable to doubly periodic parametric domains
- if parametric domain is doubly periodic, use solve_poisson_manifold_periodic()
"""
function solve_poisson_manifold(domain,partition,p,degree,u,metric_func;isperiodic=ntuple(i->false,length(partition)))

  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,partition,isperiodic=isperiodic))

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

"""
solve on manifold with doubly periodic parametric domain
use lagrange mulitplers to enforce compatibility and zeromean constraints

Note, only works in 1D (i.e. circle)
- in 2D (for sphere), the error blows up due to pole signularity
"""
function solve_poisson_manifold_periodic(domain,partition,p,degree,u,metric_func;name::ManifoldName)

  d = length(partition)
  if d > 1
    @assert name==cubedsphere "Must be 1D!"
  end

  model = CartesianDiscreteModel(domain, partition, isperiodic=ntuple(x->true,d))
  Ω = Triangulation(model)
  m = Metric(metric_func,Ω)

  if name == cubedsphere
    println("cubed sphere!")
    m = Metric(cubedsphere,Ω)
  end

  dΩ = Measure(Ω,degree)
  dΩg =  Measure(m,Ω,degree)

  # check zero mean
  println("Zero mean: ", sum(∫(u)dΩ) )

  # check compatibility
  ucf = CellField(u,Ω)
  println("Compatibility: ", sum(∫( surface_laplacian(ucf,m))dΩ) )



  ################################################################################
  #### Method 4
  #### Mixed form -- with lagrange multiplers
  ################################################################################
  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:L2)
  U = TrialFESpace(V)

  T = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p); conformity=:Hdiv)
  S = TrialFESpace(T)

  Λ = ConstantFESpace(model)
  M = TrialFESpace(Λ)

  X = MultiFieldFESpace([S,U,M])
  Y = MultiFieldFESpace([T,V,Λ])

  ### Sigma_exact
  sigma_ex(x) = surface_gradient(u,m)(x)
  _X = MultiFieldFESpace([S,M])
  _Y = MultiFieldFESpace([T,Λ])

  biformS((s,μ),(t,λ)) = ∫( s⋅t )dΩg + ∫( surface_divergence(s,m)*λ )dΩg  + ∫( surface_divergence(t,m)*μ )dΩg
  liformS((t,λ)) = ∫( sigma_ex ⋅ t )dΩg
  op = AffineFEOperator(biformS,liformS,_X,_Y)
  sigma_exh,μh = solve(LUSolver(),op)

  ### dual form
  _rhs = -1.0*surface_divergence(sigma_exh,m)

  biformX((s,u,μ),(t,v,λ)) = (  ∫( s⋅t  )dΩg + ∫( wave_divergence(t,m)*u )dΩ
                              + ∫( surface_divergence(s,m)*v  )dΩg
                              + ∫(v*μ)dΩg + ∫(λ*u)dΩg
                        )
  liformY((t,v,λ)) = ∫( -(_rhs*v) )dΩg  + ∫(λ*u)dΩg

  op = AffineFEOperator(biformX,liformY,X,Y)
  sh,uh,μh = solve(LUSolver(),op)

  #### Compute errors
  e = sum(∫((uh-u)⊙(uh-u))dΩ)
  eg = sum(∫((uh-u)⊙(uh-u))dΩg)

  return e, eg

end


"""
solve wave equation on manifold with doubly periodic parametric space

Solve
   ∫ u⋅v dΩg - ∫ p* wave_div() dΩ = ∫ vfᵤ dΩg
  ∫ pq dΩg + ∫ q* surf_div(u) dΩg = ∫ vfₚ dΩg
"""
function solve_wave_manifold(manifold,n,p,degree,uex,pex)
  _metric_func = metrics[manifold]
  domain = domains[manifold]

  d = Int(length(domain)/2)

  model = CartesianDiscreteModel(domain, ntuple(x->n,d),isperiodic=ntuple(x->true,d))

  Ω = Triangulation(model)
  m = Metric(_metric_func,Ω)


  dΩ = Measure(Ω, degree)
  dΩg =  Measure(m,Ω,degree)

  ucf = CellField(uex,Ω)
  pcf = CellField(pex,Ω)

  u0 = ucf + surface_gradient(pcf,m)
  p0 = pcf + surface_divergence(ucf,m)

  ### FE problem
  V = TestFESpace(model,ReferenceFE(raviart_thomas,Float64,p),conformity=:HDiv)
  U = TrialFESpace(V)

  Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,p),conformity=:L2)
  P = TrialFESpace(Q)

  X = MultiFieldFESpace([U,P])
  Y = MultiFieldFESpace([V,Q])

  # Weak formulation
  wave_biform((u,p),(v,q)) = ( ∫( u⋅v )dΩg
                    + ∫( -1.0*p*(wave_divergence(v,m))   )dΩ
                    + ∫( p*q )dΩg
                    + ∫( (surface_divergence(u,m))*q )dΩg
                    )
  wave_liform((v,q)) = ∫( u0⋅v + p0*q  )dΩg

  op = AffineFEOperator(wave_biform,wave_liform,X,Y)
  uh, ph = solve(LUSolver(),op)

  # Error
  eu =  sum(∫((uh-uex)⊙(uh-uex))dΩ)
  ep = sum(∫((ph-pex)⊙(ph-pex))dΩ)

  eu_g =  sum(∫((uh-uex)⊙(uh-uex))dΩg)
  ep_g = sum(∫((ph-pex)⊙(ph-pex))dΩg)

  return eu,ep,eu_g,ep_g
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
