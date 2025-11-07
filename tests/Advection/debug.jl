using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed

using DrWatson


using Gridap.CellData, Gridap.Geometry

include("../convergence_tools.jl")
include("advection_funcs.jl")
vX = panel_to_cartesian(tangent_vec(vecX))
u = panel_to_cartesian(u0)
uvX = panel_to_cartesian(u0vecX)


function my_mean( Bu_n::SkeletonPair)
plus  = ( Bu_n.plus)
minus = ( Bu_n.minus)
0.5*( plus - minus  )
end

function _my_mean(j::SkeletonPair,vel::CellField,u::CellField)
0.5*( (j.plus⋅vel.plus)*u.plus + (j.minus⋅vel.minus)*u.minus )
end

function _my_other_mean(j::SkeletonPair,vel::CellField,u::CellField,meas)
0.5*( (j.plus⋅vel.plus)*u.plus*meas.plus + (j.minus⋅vel.minus)*u.minus*meas.minus )
end

function my_jump(j::SkeletonPair,ginv::SkeletonPair,n::SkeletonPair,u::CellField)
u.plus*(j.plus⋅(ginv.plus⋅n.plus) ) + u.minus*(j.minus⋅(ginv.minus⋅n.minus) )
end
############ debug
# nprocs = 2
# ranks = with_debug() do distribute
#   distribute(LinearIndices((nprocs,)))
# end


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))
models = get_octree_refined_models(ranks,3)
# models = get_refined_models(2)
# models =  get_distributed_refined_models(ranks,nprocs,2)

panel_model = models[1]

das = FullyAssembledRows()

p_fe = 1
panel_ids = get_panel_ids(panel_model)
degree = 2*(p_fe + 1)

Ω_panel = Triangulation(das,panel_model)
dΩ = Measure(Ω_panel,degree)

# Λ = SkeletonTriangulation(panel_model)
Λ = SkeletonTriangulation(das,panel_model)
dΛ = Measure(Λ,degree)
n_Λ = get_normal_vector(Λ)


# using Gridap.Helpers, Gridap.Geometry, Gridap.CellData, Gridap.Fields
# function GridapGeosciences.panelwise_cellfield(f::Function,trian::BodyFittedTriangulation,panel_ids::AbstractArray{Int})
#   println("new cell fields")
#   @check length(panel_ids) == num_cells(trian) "\n Incorrect panel ids"
#   cell_field = map(p->GenericField(f(p)),panel_ids)
#   CellData.GenericCellField(cell_field,trian,PhysicalDomain())

#   _cf = _cell_data(f,trian,panel_ids)
#   CellData.GenericCellField(_cf,trian,ReferenceDomain())

# end

# function _cell_data(f::Function,trian,panel_ids::AbstractArray)
#   # make physical cf
#   cell_field = map(p->GenericField(f(p)),panel_ids)
#   cf = CellData.GenericCellField(cell_field,trian,PhysicalDomain())

#   panel_model = get_background_model(trian)
#   cmap = get_cell_map(get_grid(panel_model))

#   # compose with the map from refernce cell -> physical cell
#   lazy_map(∘,get_data(cf),cmap)

# end

_rhs(p) = αβ -> u(p)(αβ) + surfdiv(contra_v(uvX))(p)(αβ)

v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)
rhs_cf = panelwise_cellfield(_rhs,Ω_panel,panel_ids)

Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)

# hard code RT space as order 1 -- for velocity
V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
U = TrialFESpace(V)

vel = v_contr_cf

meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
meas_cf_skel = panelwise_cellfield(sqrtg,Λ)
_meas_cf = CellField(_sqrtg,Ω_panel)

a_Ω(u,v) = ∫( (u*v)*meas_cf )dΩ - ∫( (u*(∇(v)⋅vel) )*meas_cf )dΩ

### volume stabilisation term
a_s1(u,v) = ∫( my_mean((vel*u)⋅n_Λ)*jump(v)*meas_cf_skel.plus   )dΛ
# a_s1(u,v) = ∫( my_mean((vel*u)⋅n_Λ)*jump(v)*meas_cf   )dΛ

# jac_cf = panelwise_cellfield(forward_jacobian,Λ)
# ginv_cf = panelwise_cellfield(inv_metric,Λ)
# a_s1(u,v) = ∫( _my_mean(jac_cf,vel,u)⋅my_jump(jac_cf,ginv_cf,n_Λ,v)*meas_cf_skel.plus   )dΛ
# a_s1(u,v) = ∫( _my_other_mean(jac_cf,vel,u,meas_cf_skel)⋅my_jump(jac_cf,ginv_cf,n_Λ,v)   )dΛ

### upwinding stabilisation term
upwind = abs( (vel⋅ n_Λ).plus)
a_s2(u,v) = ∫(  0.5*(upwind)*jump(u)*jump(v)*meas_cf_skel.plus  )dΛ
# a_s2(u,v) = ∫(  0.5*(upwind)*jump(u)*jump(v)*meas_cf   )dΛ

# cell_geo_map = geo_map_func(panel_ids)
# n = pushforward_normal(Λ,cell_geo_map)
# n = pushforward_normal(Λ)
# a_s2(u,v) = ∫(  (0.5*(upwind)*jump(u*n)⋅jump(v*n))*meas_cf_skel.plus   )dΛ


biform_advection(p,q) =  a_Ω(p,q) + a_s1(p,q) + a_s2(p,q)
liform_advection(q) = ∫( (rhs_cf*q)*meas_cf )dΩ

assem = SparseMatrixAssembler(P,Q,das)
op = AffineFEOperator(biform_advection,liform_advection,P,Q,assem)

ls = LUSolver()
# uh = solve(ls,op)
A = get_matrix(op)
b = get_vector(op)
ns = numerical_setup(symbolic_setup(ls,A),A)
x = Gridap.Algebra.allocate_in_domain(A); fill!(x,0.0)
solve!(x,ns,b)
uh = FEFunction(P,x)


eu = l2((uh-u_cf)*meas_cf,dΩ)

i_am_main(ranks) && println(num_cells(panel_model))
i_am_main(ranks) && println(eu)

dir = datadir("Advection_test")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)
# cell_geo_map = geo_map_func(Ω_panel)
cell_geo_map = geo_map_func(get_panel_ids(Ω_panel))
labels = ["uh","u","eu"]
panel_cfs = [uh,u_cf,uh-u_cf]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/advection", cellfields=cellfields,append=false,geo_map=cell_geo_map)
