
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

_W = TestFESpace(model,
  ReferenceFE(lagrangian, Float64, p),
  conformity=:L2)
R = TrialFESpace(_W)

V = TestFESpace(model,
                ReferenceFE(raviart_thomas,Float64,p),
                conformity=:Hdiv)
U = TrialFESpace(V)

β_cf = CellField(β,Λ)
pts = get_cell_points(Λ)
(β_cf.plus ⋅ n_Λ.plus)(pts)



# project velocity onto Hdiv
_F = interpolate(β,U)
rhs_cf = CellField(rhs,Ω)

writevtk(Λ,dir*"/flat_test.vtu", cellfields=["plus"=>abs( (_F⋅n_Λ).plus),"minus"=>abs((_F⋅n_Λ).minus)],append=false)


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
# panel_model = Gridap.Adaptivity.refine(panel_model)

panel_ids = get_panel_ids(panel_model)

p_fe = 1
degree = 2*(p_fe + 1)
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,degree)

Λ = SkeletonTriangulation(panel_model)
dΛ = Measure(Λ,degree)
n_Λ = get_normal_vector(Λ)
pts = get_cell_points(Λ)

_sqrtg_cf = CellField(sqrtg,Λ)
sqrtg_cf = change_domain(_sqrtg_cf,PhysicalDomain(),ReferenceDomain())


panel_cfs = [sqrtg_cf.plus, sqrtg_cf.minus, sqrtg_cf.minus-sqrtg_cf.plus]
labels = ["g_plus", "g_minus", "diff"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)

skel_panel_ids = get_panel_ids(Λ)
skel_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skel_panel_ids.plus)
writevtk(Λ,dir*"/advection",cellfields=cellfields,append=false,geo_map=skel_geo_map)






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
writevtk(Ω_panel,dir*"/advection", cellfields=cellfields,append=false,geo_map=cell_geo_map)


Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)

V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
U = TrialFESpace(V)

# vel = interpolate(v_contr_cf,U)

_vel = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
vel = interpolate(_vel,U)
upwind = abs((vel⋅ n_Λ).plus)
upwind(pts)

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
labels = ["upwind_plus","upwind_minus","upwind_diff"]
panel_cfs = [abs((vel⋅ n_Λ).plus),abs((vel⋅ n_Λ).minus),abs((vel⋅ n_Λ).minus)-abs((vel⋅ n_Λ).plus)]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Λ,dir*"/advection", cellfields=cellfields,append=false,geo_map=skel_geo_map)




# v_skel = panelwise_cellfield(contra_v(vX),Λ.trian.plus,skel_panel_ids.plus)
# v_skel = interpolate(_v_skel.plus,U)
# _upwind = v_skel.plus ⋅ n_Λ.plus
# _upwind(pts)

# (upwind - _upwind)(pts)

# upwind(pts)./1

# cdata = get_data(upwind)
# _cdata = lazy_map(Operation(_my_sign),cdata)
# upwind_sign = GenericCellField(_cdata,Λ,ReferenceDomain())
# upwind_sign(pts)./1
# upwind_sign = Operation(my_sign)(upwind)
# upwind(pts)
# upwind_sign(pts)

function _my_sign(Fn)
  c = 0.5
  if Fn < eps()
    return 0.0
  elseif Fn < 0.0
    return  -0.5
  end
  return c
end


function my_upwind(a::SkeletonPair)
  Fn = a.plus
  c = Operation(my_sign)(Fn)

  return c * Fn
end

function my_mean( Bu_n::SkeletonPair, sqrtg_cf::CellField)
  plus  = ( Bu_n.plus)#*sqrtg_cf.plus
  minus = ( Bu_n.minus)#*sqrtg_cf.minus
  0.5*( plus - minus  )
end

meas_cf = CellField(sqrtg,Ω_panel)
a_Ω(u,v) = ∫( (u*v)*meas_cf )dΩ - ∫( (u*(∇(v)⋅vel) )*meas_cf )dΩ
# A1 = assemble_matrix(a_Ω,P,Q)


a_s1(u,v) = ∫( my_mean((vel*u)⋅n_Λ, sqrtg_cf)*jump(v)*meas_cf   )dΛ
# A2 = assemble_matrix(a_s1,P,Q)

upwind = abs((vel⋅ n_Λ).plus)
# a_s2(u,v) = ∫(  my_upwind(vel⋅n_Λ)*jump(u)*jump(v)*sqrtg_cf.plus   )dΛ
a_s2(u,v) = ∫(  0.5*(upwind)*jump(u)*jump(v)*meas_cf   )dΛ
# A3 = assemble_matrix(a_s2,P,Q)


a(u,v) =  a_Ω(u,v) + a_s1(u,v) + a_s2(u,v)
l(v) = ∫( (rhs_cf*v)*meas_cf )dΩ

op = AffineFEOperator(a,l,P,Q)
ph = solve(LUSolver(),op)

l2((ph-p_cf)*meas_cf,dΩ)

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
labels = ["ph","p0","ep"]
panel_cfs = [ph,p_cf,ph-p_cf]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/advection", cellfields=cellfields,append=false,geo_map=cell_geo_map)
