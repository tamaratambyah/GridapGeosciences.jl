using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Test
using GridapPETSc

dir = datadir("Hcurl")
!isdir(dir) && mkdir(dir)

include("../../convergence_tools.jl")
include("../../Geophysical/Williamson2Test.jl")
include(srcdir("Helpers/overloads.jl"))
# include("../../Geophysical/CurlConformingFESpacesFixes.jl")

## pullback 3D vector to 3D chart
inv_jacobian(p) = x -> inv(forward_jacobian_3D(p)(x))
contra_v_3D(vecX::Function,p::Int) = αβ -> inv_jacobian(p)(αβ) ⋅ vecX(p)(αβ)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


#################### sphere
o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
        num_horizontal_uniform_refinements=0,
        num_vertical_uniform_refinements=0)
panel_model = o3model.parametric_dmodel

tags = ["top_boundary", "bottom_boundary"]

panel_ids = get_panel_ids(panel_model)
Ω = Triangulation(panel_model)
dΩ = Measure(Ω,6)
dΩ_error = Measure(Ω,8)


# ## normal vector in the chart
# function fV(p)
#   function f(γαβ)
#     xyz = forward_map_3D(p)(γαβ)
#     VectorValue(xyz[1],xyz[2],xyz[3])
#   end
# end

#### Consider a more complicated vector field
function fV(p)
  function f(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    VectorValue(0.0,xyz[3],xyz[1]^2)
  end
end
vec_contra_cf = panelwise_cellfield(contra_v_3D(fV),Ω,panel_ids)

### lowest order nedelec
order = 0
value_type = Float64
reffe =  ReferenceFE(nedelec,Float64,order)
R = TestFESpace(panel_model,reffe;conformity=:Hcurl,dirichlet_tags=tags)
H = TrialFESpace(R,vec_contra_cf)

nedelec_reffe = Gridap.ReferenceFEs.NedelecRefFE(value_type,HEX,order)


ndofs = Gridap.ReferenceFEs.num_dofs(nedelec_reffe)
face_own_dofs  = Gridap.ReferenceFEs.get_face_own_dofs(nedelec_reffe)
println(face_own_dofs)

########### Interpolation
# vec_contra_h = interpolate(vec_contra_cf,H)
# d = collect(get_cell_dof_values(vec_contra_h.fields.item_ref[]))
# _d = map(x->round.(x;digits=8),d)

############ Evaluation on dof basis
##### by evaluating on the dof basis, we avoid noise of the scatter/gather in the
##### FE space
dof_basis = get_fe_dof_basis(H)
vec_contra_h = dof_basis(vec_contra_cf)
d = collect(vec_contra_h.item_ref[])
_d = map(x->round.(x;digits=8),d)

##### Panel 6 is the master, panel 1 is the slave
## Edge [1 2]:  Panel 6 dof 1 -> Panel 1 dof 1
_d[6][1]
_d[1][1]

## Edge [5 6]:  Panel 6 dof 2 -> Panel 1 dof 3
_d[6][2]
_d[1][3]

## Edge [1 5]:  Panel 6 dof 5 -> Panel 1 dof 9
_d[6][5]
_d[1][9]

## Edge [2 6]:  Panel 6 dof 7 -> Panel 1 dof 11
_d[6][7]
_d[1][11]

####### Panel 3 is the master, panel 1 is the slave
## Edge [3 4]:  Panel 3 dof 1 -> Panel 1 dof 2
_d[3][1]
_d[1][2]

## Edge [7 8]:  Panel 3 dof 3 -> Panel 1 dof 4
_d[3][3]
_d[1][4]

#### WRONG!
## Edge [3 7]:  Panel 3 dof 9 -> Panel 1 dof 10
_d[3][9]
_d[1][10]

#### WRONG!
## Edge [4 8]:  Panel 3 dof 11 -> Panel 1 dof 12
_d[3][11]
_d[1][12]

metric_cf = panelwise_cellfield(metric,Ω,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω,panel_ids)
covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω,panel_ids)

_e = vec_contra_cf - vec_contra_h
el2_interp =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ_error))



############ L2 projection
a(u,v) = ∫( u⋅( metric_cf⋅v)*meas_cf )dΩ
l(v) = ∫( vec_contra_cf⋅(metric_cf⋅v)*meas_cf )dΩ
op = AffineFEOperator(a,l,H,R)
vec_contra_h = solve(LUSolver(),op)

d = collect(get_cell_dof_values(vec_contra_h.fields.item_ref[]))
_d = map(x->round.(x;digits=8),d)

##### Panel 6 is the master, panel 1 is the slave
## Edge [1 2]:  Panel 6 dof 1 -> Panel 1 dof 1
_d[6][1]
_d[1][1]

## Edge [5 6]:  Panel 6 dof 2 -> Panel 1 dof 3
_d[6][2]
_d[1][3]

## Edge [1 5]:  Panel 6 dof 5 -> Panel 1 dof 9
_d[6][5]
_d[1][9]

## Edge [2 6]:  Panel 6 dof 7 -> Panel 1 dof 11
_d[6][7]
_d[1][11]

####### Panel 3 is the master, panel 1 is the slave
## Edge [3 4]:  Panel 3 dof 1 -> Panel 1 dof 2
_d[3][1]
_d[1][2]

## Edge [7 8]:  Panel 3 dof 3 -> Panel 1 dof 4
_d[3][3]
_d[1][4]

### WRONG!!!!
## Edge [3 7]:  Panel 3 dof 9 -> Panel 1 dof 10
_d[3][9]
_d[1][10]

#### WRONG
## Edge [4 8]:  Panel 3 dof 11 -> Panel 1 dof 12
_d[3][11]
_d[1][12]


vec_proj_h = covariant_basis_cf ⋅vec_contra_h

_e = vec_contra_cf - vec_contra_h
el2_proj =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ))
