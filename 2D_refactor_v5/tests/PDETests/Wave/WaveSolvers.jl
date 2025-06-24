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
