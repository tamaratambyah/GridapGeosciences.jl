"""
solve helmholtz equation on doubly periodic, flat domain/
Consider uex that fullfil compatibility and zeromean constraints
"""
function helmholtz_periodic(domain,partition,p,degree,uex)

  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,partition,isperiodic=ntuple(x->true,length(partition))))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  ucf = CellField(uex,Ω)
  compat = ( sum( ∫( ucf)dΩ) + sum(∫( laplacian(ucf))dΩ  ) )
  println("Compatibility: ", compat)

  rhs = (ucf + laplacian(ucf))

    ### FE problem -- single field
  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
  U = TrialFESpace(V)


  poisson_biform(u,v) = ∫(u*v)dΩ -  ∫( gradient(u)⋅gradient(v)  )dΩ
  poisson_liform(v) = ∫(  rhs*v )dΩ

  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
  uh = solve(LUSolver(),op)

  return sum(∫((uh-uex)⊙(uh-uex))dΩ)
end

"""
solve helmholtz equation on doubly periodic, flat domain, in dual form
Consider uex that fullfil compatibility and zeromean constraints
"""
function helmholtz_dual_periodic(domain,partition,p,degree,uex)

  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,partition,isperiodic=ntuple(x->true,length(partition))))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  ucf = CellField(uex,Ω)
  compat = ( sum( ∫( ucf)dΩ) + sum(∫( laplacian(ucf))dΩ  ) )
  println("Compatibility: ", compat)

  rhs = (ucf + laplacian(ucf))

  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:L2)
  U = TrialFESpace(V)

  T = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p); conformity=:Hdiv)
  S = TrialFESpace(T)

  X = MultiFieldFESpace([S,U])
  Y = MultiFieldFESpace([T,V])

  biformX((s,u),(t,v)) = (  ∫( s⋅t)dΩ + ∫( divergence(t)*u )dΩ
                          + ∫( u*v )dΩ   + ∫( divergence(s)*v  )dΩ
                        )
  liformY((t,v)) = ∫( rhs*v )dΩ

  op = AffineFEOperator(biformX,liformY,X,Y)
  sh,uh = solve(LUSolver(),op)


  return sum(∫((uh-uex)⊙(uh-uex))dΩ)
end



"""
solve helmholtz equation on doubly periodic, parametric domain  domain.
Consider uex that fullfil compatibility and zeromean constraints
"""
function helmholtz_manifold_periodic(manifold,n,p,degree,uex)
  _metric_func = metrics[manifold]
  domain = domains[manifold]

  d = Int(length(domain)/2)

  model = CartesianDiscreteModel(domain, ntuple(x->n,d), isperiodic=ntuple(x->true,d))
  Ω = Triangulation(model)
  m = Metric(_metric_func,Ω)

  dΩ = Measure(Ω,degree)
  dΩg =  Measure(m,Ω,degree)


  ucf = CellField(uex,Ω)
  # check compatibility
  compat = sum( ∫(ucf)dΩ + ∫( surface_laplacian(ucf,m))dΩ  )
  println("Compatibility: ", compat)

  rhs = ucf + 1.0*(surface_laplacian(ucf,m))

  #### FE Problem -- no lagrange multiplers
  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
  U = TrialFESpace(V)

  poisson_biform(u,v) = ∫(u*v)dΩg -  ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
  poisson_liform(v) = ∫(  rhs*v )dΩg
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

  uh = solve(LUSolver(),op)

  #### Compute errors
  e = l2(uh-ucf,dΩ)
  eg = l2(uh-ucf,dΩg)
  return e,eg

end
