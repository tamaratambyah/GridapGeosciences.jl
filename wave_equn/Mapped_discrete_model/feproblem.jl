
function fe_problem(model,p,u_ex)
  degree = 2*(p+1)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  lg_reffe = ReferenceFE(lagrangian, Float64, p)
  V = FESpace(model, lg_reffe; conformity=:L2,dirichlet_tags="boundary")
  U = TrialFESpace(V,u_ex)

  a(u,v) = ∫( u⊙v )dΩ
  b(v) = ∫( u_ex*v )dΩ

  op = AffineFEOperator(a,b,U,V)

  return uh = solve(LUSolver(),op)
end



function fe_problem_Hdiv(model,p,u_ex)
  degree = 2*(p+1)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  rt_reffe = ReferenceFE(raviart_thomas, Float64, p)
  V = FESpace(model, rt_reffe; conformity=:Hdiv,dirichlet_tags="boundary")
  U = TrialFESpace(V,u_ex)

  a(u,v) = ∫( u⋅v )dΩ
  b(v) = ∫( u_ex⋅v )dΩ

  op = AffineFEOperator(a,b,U,V)

  return uh = solve(LUSolver(),op)
end
