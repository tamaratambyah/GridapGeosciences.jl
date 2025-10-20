using Gridap
using GridapGeosciences
using DrWatson
using FillArrays
using Gridap.Geometry, Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Test


using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using MPIPreferences
MPIPreferences.use_jll_binary()


dir = datadir("Distributed")
include("../convergence_tools.jl")

s_models = get_refined_models(3)

model = s_models[1]
Ω_panel = Triangulation(model)
panel_ids = get_panel_ids(Ω_panel)


cell_geo_map = latlon_geo_map_func(get_panel_ids(Ω_panel))
writevtk(Ω_panel,dir*"/ambient_model", append=false,ascii=true,geo_map=cell_geo_map)

function latlon_geo_map_func(panel_ids::AbstractArray{Int})
  println("latolon serial geo map")
  return lazy_map(p -> Cartesian2SphereicalMap() ∘ MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
end





################################################################################
#### Octree
################################################################################
MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

parametric_octree_dmodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=3)

panel_model = parametric_octree_dmodel.parametric_dmodel

trian = Triangulation(panel_model)
get_panel_ids(panel_model)
get_owned_panel_ids(panel_model)

n_ref_lvls = 3
coarse_model = true

dmodels = get_octree_refined_models(ranks,n_ref_lvls,coarse_model)

map(x->num_cells(x),dmodels)

s_models = get_refined_models(n_ref_lvls,coarse_model)
map(x->num_cells(x),s_models)

map(x->get_panel_ids(x),dmodels)


################################################################################
#### Advection dg debugging
################################################################################
include("../Advection/advection_funcs.jl")
vX = panel_to_cartesian(tangent_vec(vecX))
u = panel_to_cartesian(u0)
uvX = panel_to_cartesian(u0vecX)

# panel_model = s_models[1]
panel_model = dmodels[2]
p_fe = 1

panel_ids = get_panel_ids(panel_model)
degree = 2*(p_fe + 1)

Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,degree)



# function Gridap.Geometry.SkeletonTriangulation(model::DistributedParametricDiscreteModel;kwargs...)
#   println("my skeleton -- default is with_ghost")
#   SkeletonTriangulation(with_ghost,model;kwargs...)
# end

# Λ = SkeletonTriangulation(with_ghost,panel_model)
Λ = SkeletonTriangulation(panel_model)
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

## albertos tests
uh = interpolate(u_cf,P)
cf = jump(uh)

dc = ∫( jump(uh) )dΛ
sum(dc)

vel = interpolate(v_contr_cf,U)
meas_cf = CellField(sqrtg,Ω_panel)

a_Ω(u,v) = ∫( (u*v)*meas_cf )dΩ - ∫( (u*(∇(v)⋅vel) )*meas_cf )dΩ

# function my_mean( Bu_n::SkeletonPair)
#   plus  = ( Bu_n.plus)
#   minus = ( Bu_n.minus)
#   0.5*( plus - minus  )
# end

function _my_mean(j::SkeletonPair,vel::CellField,u::CellField)
  0.5*( (j.plus⋅vel.plus)*u.plus + (j.minus⋅vel.minus)*u.minus )
end

function my_jump(j::SkeletonPair,ginv::SkeletonPair,n::SkeletonPair,u::CellField)
  u.plus*(j.plus⋅(ginv.plus⋅n.plus) ) + u.minus*(j.minus⋅(ginv.minus⋅n.minus) )
end

### volume stabilisation term
# a_s1(u,v) = ∫( my_mean((vel*u)⋅n_Λ)*jump(v)*meas_cf   )dΛ

jac_cf = panelwise_cellfield(forward_jacobian,Λ)
ginv_cf = panelwise_cellfield(_analytic_inv_metric,Λ)
a_s1(u,v) = ∫( _my_mean(jac_cf,vel,u)⋅my_jump(jac_cf,ginv_cf,n_Λ,v)*meas_cf   )dΛ

upwind = abs((vel⋅ n_Λ).plus)
### upwinding stabilisation term
# a_s2(u,v) = ∫(  0.5*abs((vel⋅ n_Λ).plus)*jump(u)*jump(v)*meas_cf   )dΛ

cell_geo_map = geo_map_func(panel_ids)
n = pushforward_normal(Λ,cell_geo_map)


function GridapGeosciences.pushforward_normal(trian::Gridap.Geometry.TriangulationView)
  # cf = _pushforward_normal(trian.parent)
#
  # data = get_data(cf)
  # _data = Gridap.Geometry.restrict(data, trian.cell_to_parent_cell)
  # CellData.GenericCellField(_data,trian,DomainStyle(cf))
   _pushforward_normal(trian.parent)
end

n = pushforward_normal(Λ)


a_s2(u,v) = ∫(  0.5*(upwind)*jump(u*n)⋅jump(v*n)*meas_cf   )dΛ

biform_advection(p,q) =  a_Ω(p,q) + a_s1(p,q) + a_s2(p,q)
liform_advection(q) = ∫( (rhs_cf*q)*meas_cf )dΩ

op = AffineFEOperator(biform_advection,liform_advection,P,Q)

# uh = solve(ls,op)
ls = LUSolver()
A = get_matrix(op)
b = get_vector(op)
ns = numerical_setup(symbolic_setup(ls,A),A)
x = Gridap.Algebra.allocate_in_domain(A); fill!(x,0.0)
solve!(x,ns,b)
uh = FEFunction(P,x)


eu = l2((uh-u_cf)*meas_cf,dΩ)

lvl = nref(nc(panel_model))
cell_geo_map = geo_map_func(Ω_panel)
labels = ["uh","u","eu"]
panel_cfs = [jump(uh),u_cf,uh-u_cf]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe", cellfields=cellfields,append=false,geo_map=cell_geo_map)
