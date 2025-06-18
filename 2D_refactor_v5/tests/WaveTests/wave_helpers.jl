
function solve_wave_manifold(domain,partition,p,degree,metric_func,uex,pex;name::ManifoldName)
  model = CartesianDiscreteModel(domain, partition,isperiodic=ntuple(x->true,length(partition)))

  Ω = Triangulation(model)
  m = Metric(metric_func,Ω)
  if name == cubedsphere
    println("cubed sphere!")
    m = Metric(cubedsphere,Ω)
  end

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
