
################################################################################
### FLAT TEST
################################################################################
using Gridap
β(x) = VectorValue(1.0,0.0)
u0(x) = 1 + cos(2*π*x[1])
uβ(x) = VectorValue(1 + cos(2*π*x[1]) ,0.0)

rhs(x) = u0(x) + divergence(uβ)(x)


CFL = 0.1
n = 16
p = 1
degree = 3*(p+1)

L = 1.0
domain = (0.0, L, 0.0, L)
partition = (n, n)
model = CartesianDiscreteModel(domain, partition, isperiodic=(true,true))

Ω = Triangulation(model)
dΩ = Measure(Ω, degree)
Λ = SkeletonTriangulation(model)
dΛ = Measure(Λ,degree)
n_Λ = get_normal_vector(Λ)

W = TestFESpace(model,
  ReferenceFE(lagrangian, Float64, p),
  conformity=:L2)
R = TrialFESpace(W)

V = TestFESpace(model,
                ReferenceFE(raviart_thomas,Float64,p),
                conformity=:Hdiv)
U = TrialFESpace(V)

# project velocity onto Hdiv
F = interpolate(β,U)
rhs_cf = CellField(rhs,Ω)

function my_sign(Fn)
  c = 0.5
  if Fn[1] < 0.0
    c = -0.5
  end

  return c
end

a_Ω(u,v) = ∫( u*v )dΩ - ∫( u*(∇(v)⋅F)  )dΩ
a_s(u,v) = ∫( mean(F*u)⋅jump(v*n_Λ ) +  (my_sign∘( (F⋅n_Λ).plus )*( (F⋅n_Λ).plus ) )*jump(u*n_Λ )⋅jump(v*n_Λ )   )dΛ
a(u,v) =  a_Ω(u,v) + a_s(u,v)

l(v) = ∫( rhs_cf*v )dΩ

op = AffineFEOperator(a,l,R,W)
uh = solve(LUSolver(),op)

l2(w,dΩ) = sum(∫(  w⊙w)dΩ)
l2((uh-u0),dΩ)

writevtk(Ω,dir*"/flat_test" * ".vtu", cellfields=["uh"=>uh,"u"=>u0,"eu"=>uh-u0,"rhs"=>rhs],append=false)




################################################################################
### MANIFOLD TEST
################################################################################
panel_model = coarse_parametric_model()
panel_model = Gridap.Adaptivity.refine(panel_model)
panel_model = Gridap.Adaptivity.refine(panel_model)

p_fe = 1
degree = 2*(p_fe + 1)
Ω_panel = Triangulation(panel_model)
panel_ids = get_panel_ids(panel_model)
dΩ = Measure(Ω_panel,degree)

Λ = SkeletonTriangulation(panel_model)
dΛ = Measure(Λ,degree)
n_Λ = get_normal_vector(Λ)


vecX(XYZ) = VectorValue(-XYZ[2],XYZ[3],0.0)
h0(XYZ) = exp(-(XYZ[2]^2 + XYZ[3]^2)  )
h0v(XYZ) = VectorValue(-XYZ[2]*exp(-(XYZ[2]^2 + XYZ[3]^2)  ),  XYZ[3]*exp(-(XYZ[2]^2 + XYZ[3]^2)  ),  0.0)


vX = panel_to_cartesian(tangent_vec(vecX))
φ = panel_to_cartesian(h0)
_vX = panel_to_cartesian(h0v)
_rhs(p) = αβ -> φ(p)(αβ) + surfdiv(contra_v(_vX))(p)(αβ)



v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
p_cf = panelwise_cellfield(φ,Ω_panel,panel_ids)
rhs_cf = panelwise_cellfield(_rhs,Ω_panel,panel_ids)

sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(_vX)),Ω_panel,panel_ids)
_rhs_cf = p_cf + sdiv_cf

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
labels = ["rhs","p","v","_rhs", "er"]
panel_cfs = [rhs_cf,p_cf,v_contr_cf,_rhs_cf, rhs_cf-_rhs_cf]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/solT_0" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)


Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)

V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
U = TrialFESpace(V)


# vel = interpolate(v_contr_cf,U)
vel = v_contr_cf
meas_cf = CellField(sqrtg,Ω_panel)
_meas_cf = CellField(sqrtg,Λ)
function my_sign(Fn)
  c = 0.5
  if Fn[1] < 0.0
    c = -0.5
  end
  return c
end

a_Ω(u,v) = ∫( (u*v)*meas_cf )dΩ - ∫( (u*(∇(v)⋅vel) )*meas_cf )dΩ
a_s(u,v) = ∫( (mean(vel*u)⋅jump(v*n_Λ ))*_meas_cf +  ((my_sign∘( ((vel)⋅n_Λ).plus )*( ((vel)⋅n_Λ).plus ) )*jump(u*n_Λ )⋅jump(v*n_Λ ))*_meas_cf   )dΛ
a(u,v) =  a_Ω(u,v) + a_s(u,v)

l(v) = ∫( (rhs_cf*v)*meas_cf )dΩ

op = AffineFEOperator(a,l,P,Q)
ph = solve(LUSolver(),op)

l2((ph-p_cf)*meas_cf,dΩ)

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
labels = ["ph","p0","ep"]
panel_cfs = [ph,p_cf,ph-p_cf]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/solT_0" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)




### Transient

_mass(t,p,q) = mass(p,q)
_res(t,p,q) =  res(p,q)
jac(t,p,dp,q) = res(dp,q)
jac_t(t,p,dtp,q) = mass(dtp,q)
# opT = TransientSemilinearFEOperator(_mass, _res, P, Q, constant_mass=true)
opT = TransientSemilinearFEOperator(_mass,_res, (jac, jac_t), P, Q, constant_mass=true)

t0, tF = 0.0, 2*π
CFL = 0.1
_dt = dx(nc(panel_model))*CFL/p_fe
dt = 0.001

# solve with SSP RK 3
solver = RungeKutta(LUSolver(), LUSolver(), dt, :EXRK_SSP_3_3)
solT = solve(solver, opT, t0, tF, ph0)

## iterate solution
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
labels = ["ph","v"]
panel_cfs = [ph0,vel]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/solT_0" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)



for (t, ph) in solT
  println(t)
  panel_cfs = [ph]
  cellfields = map((x,y) -> x=>y, labels,panel_cfs)
  writevtk(Ω_panel,dir*"/solT_$t" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
end


make_pvd(dir,"solT",1)
