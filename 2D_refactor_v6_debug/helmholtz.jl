### helmholtz
dΩ = Measure(Ω_panel,10)


meas_cf = CellField(sqrtg,Ω_panel)
inv_metric_cf = CellField(analytic_inv_metric,Ω_panel)

V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,2); conformity=:H1)
U = TrialFESpace(V)


sum(∫( rhs_cf*meas_cf  )dΩ)


poisson_biform(u,v) = ∫(u*v*meas_cf)dΩ -  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
poisson_liform(v) = ∫(  (rhs_cf*v)*meas_cf )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(e,dΩ) = sum(∫( e⋅e )dΩ)
e =  l2(uh-f_cf,dΩ)


_uh = change_domain(uh,ReferenceDomain(),PhysicalDomain())
cf_mapped = lazy_map(Broadcasting(∘),get_data(_uh),inv_f)
uh_ambient = CellData.GenericCellField(cf_mapped,Ω_sphere,PhysicalDomain() ) # am

writevtk(Triangulation(ambient_model),dir*"/ambient_model_slap",
          cellfields=["slap"=>slap_ambient,"f"=>f_ambient, "rhs"=>rhs_ambient,
          "uh"=>uh_ambient, "e"=>uh_ambient-f_ambient
          ],append=false)
