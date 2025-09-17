include("../Lauritzen_functions.jl")

dir = datadir("Transient_advection_test_supg")
!isdir(dir) && mkdir(dir)

panel_model = coarse_parametric_model()
panel_model = Gridap.Adaptivity.refine(panel_model)

p_fe = 1
CFL = 0.1


depth = panel_to_cartesian(cosine_bell)
vX(t) = panel_to_cartesian(tangent_vec(nondivergent_velocity(t)))


panel_ids = get_panel_ids(panel_model)
degree = 2*(p_fe + 1)

Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,degree)

Λ = SkeletonTriangulation(panel_model)
dΛ = Measure(Λ,degree)
n_Λ = get_normal_vector(Λ)

Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TransientTrialFESpace(Q)


V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
U = TrialFESpace(V)


metric_cf = CellField(analytic_metric,Ω_panel)
meas_cf = CellField(sqrtg,Ω_panel)

### Initial condition
depth_cf = panelwise_cellfield(depth,Ω_panel,panel_ids)
uh0 = interpolate(depth_cf, P(0.0))

## velocity
velbiform(u,v) = ∫( u⋅(metric_cf ⋅ v)*meas_cf )dΩ
velliform(t,v) = ∫( v⋅(metric_cf⋅panelwise_cellfield(contra_v(vX(t)),Ω_panel,panel_ids) )*meas_cf )dΩ
_velliform(v) = velliform(0.0,v)
op = AffineFEOperator(velbiform,_velliform,U,V)
_Fh = solve(LUSolver(),op)

v_contr_cf =  interpolate(panelwise_cellfield(contra_v(vX(0.0)),Ω_panel,panel_ids), U)
covarient_basis_cf = panelwise_cellfield(forward_jacobian,Ω_panel,panel_ids)

vec_phys = covarient_basis_cf ⋅ v_contr_cf
Fh = covarient_basis_cf ⋅ _Fh

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
panel_cfs = [vec_phys, vec_phys-Fh, uh0 ]
labels = ["u", "eu", "depth"]

cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model",cellfields=cellfields,append=false,geo_map=cell_geo_map)

vecX(XYZ) = VectorValue(-XYZ[2],XYZ[1],0.0)
_vX = panel_to_cartesian(tangent_vec(vecX))

function get_velocity(t)
  # interpolate(panelwise_cellfield(contra_v(vX(0.0)),Ω_panel,panel_ids), U)
  interpolate(panelwise_cellfield(contra_v(_vX),Ω_panel,panel_ids), U)
end

# equation for depth: advective form
vel = interpolate(panelwise_cellfield(contra_v(_vX),Ω_panel,panel_ids), U)

a_mass_Ω(dtu,v) = ∫( (dtu*v)*meas_cf )dΩ
a_mass_s(t,dtu,v) = ∫( (dtu*(vel⋅∇(v)))*meas_cf )dΩ
a_Ω(t,u,v) = ∫( ((vel⋅∇(u))*v )*meas_cf )dΩ
a_s(t,u,v) =  ∫( ((vel⋅∇(u))*(vel⋅∇(v)) )*meas_cf )dΩ


a_mass(t,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(t,dtu,v)
res(t,u,v) =  a_Ω(t,u,v) + τ*a_s(t,u,v)
jac(t,u,du,v) = a_Ω(t,du,v) + τ*a_s(t,du,v)
jac_t(t,u,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(t,dtu,v)
opT = TransientSemilinearFEOperator(a_mass, res, (jac,jac_t), P, Q)



### flux form
# a_mass(t,dtu,v) = ∫( (dtu*v)*meas_cf )dΩ
# a_Ω(t,u,v) =   ∫( -(u*(∇(v)⋅panelwise_cellfield(contra_v(vX(t)),Ω_panel,panel_ids)) )*meas_cf )dΩ
# a_s1(t,u,v) = ∫( my_mean((panelwise_cellfield(contra_v(vX(t)),Ω_panel,panel_ids)*u)⋅n_Λ)*jump(v)*meas_cf   )dΛ
# a_s2(t,u,v) = ∫(  0.5*abs((panelwise_cellfield(contra_v(vX(t)),Ω_panel,panel_ids)⋅ n_Λ).plus)*jump(u)*jump(v)*meas_cf   )dΛ

# res(t,u,v) =  a_Ω(t,u,v) + a_s1(t,u,v) + a_s2(t,u,v)
# jac(t,u,du,v) = a_Ω(t,du,v) + a_s1(t,du,v) + a_s2(t,du,v)
# jac_t(t,u,dtu,v) = ∫( (dtu*v)*meas_cf )dΩ
# opT = TransientSemilinearFEOperator(a_mass, res, (jac,jac_t), P, Q, constant_mass=true)


ls = LUSolver()


t0, tF = 0.0, 5.0
_dt = dx(nc(panel_model))*CFL/(p_fe*3.26) # Umax = 3.26 maximum zonal wind speed
dt = floor(_dt, sigdigits=1)
τ = 0.5*dt ## stabilisation

ode_solver = RungeKutta(ls,ls,dt,:EXRK_SSP_3_3)


solT  = solve(ode_solver,opT,t0,tF,uh0)
it = iterate(solT)


writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=["uh"=>uh0,"vel"=>Fh],append=false,geo_map=cell_geo_map)

Es = Float64[]
while !isnothing(it)
  data, state = it
  t, uh = data

  Fh = get_velocity(t)
  eu = l2((uh-uh0)*meas_cf,dΩ)
  println("t = ", t, "; ", eu)
  push!(Es,eu)

  writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=["uh"=>uh,"vel"=>covarient_basis_cf⋅Fh],append=false,geo_map=cell_geo_map)

  it = iterate(solT, state)
end
make_pvd(dir,"solT",1)

plot()
plot!( dt*collect(1:length(Es)), Es, lw=2,label=false)
