vX = panel_to_cartesian(tangent_vec(vecX))
u = panel_to_cartesian(u0)
uvX = panel_to_cartesian(u0vecX)

models  = get_refined_models(n_ref_lvls)

models,  = get_distributed_refined_models(ranks,nprocs,models)

panel_model = models[1]
p_fe = 1

lvl = nref(nc(panel_model))
println("nref = $lvl")

panel_ids = get_panel_ids(panel_model)
degree = 2*(p_fe + 1)

Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,degree)

Λ = SkeletonTriangulation(panel_model)

btrian = Λ.trians.items[1].minus
# panel_model = get_background_model(btrian)
panel_ids = get_panel_ids(panel_model)
Dc = num_cell_dims(panel_model)

glue = get_glue(btrian,Val(Dc))
face_2_cell = glue.tface_to_mface
face_panel_ids = panel_ids[face_2_cell]
return face_panel_ids




dΛ = Measure(Λ,degree)
n_Λ = get_normal_vector(Λ)

_rhs(p) = αβ -> u(p)(αβ) + surfdiv(contra_v(uvX))(p)(αβ)

v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)
rhs_cf = panelwise_cellfield(_rhs,Ω_panel,panel_ids)

Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)

# hard code RT space as order 1 -- for velocity
V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
U = TrialFESpace(V)


trian = Ω_panel
f = contra_v(vX)
fields = map(trian.trians,panel_ids) do trian, pids
  println(typeof(trian.trian)<:Gridap.Geometry.BodyFittedTriangulation)
  cell_field = map(p->Gridap.Fields.GenericField(f(p)),pids)
  Gridap.CellData.GenericCellField(cell_field,trian.trian,PhysicalDomain())
end
v_contr_cf = GridapDistributed.DistributedCellField(fields,trian)
vel = interpolate(v_contr_cf,U)


_a(u,v) = ∫( u⋅v )dΩ
_l(v) = ∫( v_contr_cf ⋅ v )dΩ
op = AffineFEOperator(_a,_l,U,V)
vel = solve(LUSolver(),op)

meas_cf = CellField(sqrtg,Ω_panel)

a_Ω(u,v) = ∫( (u*v)*meas_cf )dΩ - ∫( (u*(∇(v)⋅vel) )*meas_cf )dΩ

### volume stabilisation term
# a_s1(u,v) = ∫( my_mean((vel*u)⋅n_Λ)*jump(v)*meas_cf   )dΛ

jac_cf = panelwise_cellfield(forward_jacobian,Λ)
ginv_cf = panelwise_cellfield(_analytic_inv_metric,Λ)
a_s1(u,v) = ∫( _my_mean(jac_cf,vel,u)⋅my_jump(jac_cf,ginv_cf,n_Λ,v)*meas_cf   )dΛ


### upwinding stabilisation term
upwind = abs((vel⋅ n_Λ).plus)
# a_s2(u,v) = ∫(  0.5*(upwind)*jump(u)*jump(v)*meas_cf   )dΛ

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
cell_normal = get_facet_normal(Λ,cell_geo_map)
n = get_normal_vector(Λ,cell_normal)
a_s2(u,v) = ∫(  0.5*(upwind)*jump(u*n)⋅jump(v*n)*meas_cf   )dΛ


biform_advection(p,q) =  a_Ω(p,q) + a_s1(p,q) + a_s2(p,q)
liform_advection(q) = ∫( (rhs_cf*q)*meas_cf )dΩ

op = AffineFEOperator(biform_advection,liform_advection,P,Q)

# uh = solve(ls,op)
A = get_matrix(op)
b = get_vector(op)
ns = numerical_setup(symbolic_setup(ls,A),A)
x = allocate_in_domain(A); fill!(x,0.0)
solve!(x,ns,b)
uh = FEFunction(P,x)
