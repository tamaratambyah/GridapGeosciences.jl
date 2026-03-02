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
using LinearAlgebra

dir = datadir("Hcurl")
!isdir(dir) && mkdir(dir)

include("../../convergence_tools.jl")
include("../../Geophysical/Williamson2Test.jl")
include(srcdir("Helpers/overloads.jl"))
include("../../Geophysical/CurlConformingFESpacesFixes.jl")

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
_panel_model = o3model.parametric_dmodel
panel_model = _panel_model.models.item

tags = ["top_boundary", "bottom_boundary"]

panel_ids = get_panel_ids(panel_model)
Ω = Triangulation(panel_model)
dΩ = Measure(Ω,6)
dΩ_error = Measure(Ω,8)

cmap = get_cell_map(get_grid(panel_model))
coords = HEX.vertex_coords

#### Consider a complicated vector field
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

dof_basis = get_fe_dof_basis(H)
vec_contra_h = dof_basis(vec_contra_cf)
cell_dofs = collect(vec_contra_h)


########### Panel 3 is the master, panel 1 is the slave ########################
pm = 3
ps = 1

########### Edge [3 7] #######################
## In panel 3,edge [3 7] -> dof 9
## In panel 1, edge [3 7] -> dof 10
dof_3 = cell_dofs[3][9]
dof_1 = cell_dofs[1][10]


## In panel 3, Node 3 -> (γ,α,β) = (0,-π/4,-π/4)
## In panel 1, Node 3 -> (γ,α,β) = (0,π/4,-π/4)
vm = Point(0,-π/4,-π/4)
vs = Point(0,π/4,-π/4)

function slave2master(γαβ_s,p_s,p_m)
  ## 1. push slave to ambient space
  ## 2. pull to master using inverse of master map
  xyz = ForwardMap(p_s)(γαβ_s)
  inv_m = inverse_map(ForwardMap(p_m))
  γαβ_m = evaluate(inv_m,xyz)
  return γαβ_m
end

#### Node 3: given the slave parametric coord, compute the master parametric coord
γαβ_s = evaluate(cmap[1],coords[3])
γαβ_m = slave2master(γαβ_s,ps,pm)
γαβ_m ≈ vm


#### Tangent to the edge [3, 7]
tangent_m = VectorValue(0.0,0.0,1.0)

W = inv_jacobian(ps)(γαβ_s)⋅forward_jacobian(pm)(γαβ_m)
coeffs = Matrix(W)
tangent_s = W ⋅ tangent_m
tangent_s ≈ VectorValue(0,0,1)


#### Cofficient matrix in terms of slave points
function W_matrix(γαβ_s,p_s,p_m)
  inv_jacobian(ps)(γαβ_s)⋅forward_jacobian(p_m)( slave2master(γαβ_s,p_s,p_m) )
end

W_s = W_matrix(γαβ_s,ps,pm)
coeffs_s = Matrix(W_s)
coeffs_s ≈ coeffs
