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
et = Float64
reffe =  ReferenceFE(nedelec,et,order)
R = TestFESpace(panel_model,reffe;conformity=:Hcurl,dirichlet_tags=tags)
H = TrialFESpace(R,vec_contra_cf)

dof_basis = get_fe_dof_basis(H)
vec_contra_h = dof_basis(vec_contra_cf)
cell_dofs = collect(vec_contra_h)

################################################################################
########### Edge [3 7] #######################
################################################################################

########### Panel 3 is the master, panel 1 is the slave
pm = 3
ps = 1

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
γαβ_s = Point(0,π/4,-π/4)
γαβ_m = slave2master(γαβ_s,ps,pm)
γαβ_m ≈ vm


#### Tangent to the edge [3, 7]
tangent_m = VectorValue(0.0,0.0,1.0)

W = inv_jacobian(ps)(γαβ_s)⋅forward_jacobian(pm)(γαβ_m)
coeffs = Matrix(W)
tangent_s = W ⋅ tangent_m
tangent_s ≈ VectorValue(0,0,1)


#### Cofficient matrix in terms of slave points
function W_matrix(p_s,p_m)
  function _W_matrix(γαβ_s)
    W = inv_jacobian(p_s)(γαβ_s)⋅forward_jacobian(p_m)( slave2master(γαβ_s,p_s,p_m) )
    W
  end
end

W_s = W_matrix(ps,pm)(γαβ_s)
coeffs_s = Matrix(W_s)
coeffs_s ≈ coeffs


########### integral over the edge:
##### an the edge [3,7] is [0,-π/4,-π/4] -> [0,π/4,π/4]
##### we could restrict the 3D model to the edge, but as a hack, let's create a
##### 1D model that is the line [-π/4,π/4]
model_edge = UnstructuredDiscreteModel(CartesianDiscreteModel((-π/4,π/4),(1)))
Ω_edge = Triangulation(model_edge)
dΩ_edge = Measure(Ω_edge,0) ## zero order elements -> zero order quadrature for edge values


############# PANEL 3
## In panel 3,edge [3 7] -> (γ,α,β) = (0,-π/4,β) where β = [-π/4,π/4], tangent = (0,0,1)
tangent_m = VectorValue(0,0,1)
function u_dot_t_p3(β)
  γαβ = VectorValue(0.0,-π/4,β[1])
  u = inv_jacobian(pm)(γαβ) ⋅ fV(pm)(γαβ)
  u ⋅ tangent_m
end
dc_m = ∫(  u_dot_t_p3   )dΩ_edge
dc_m.dict
sum(dc_m) ≈ dof_3


############# PANEL 1
## In panel 1, edge [3 7] -> (γ,α,β) = (0,π/4,β) where β = [-π/4,π/4], tangent = (0,0,1)
dof_1
function u_dot_t_p1(β)
  γαβ_s = VectorValue(0.0,π/4,β[1])
  γαβ_m = slave2master(γαβ_s,ps,pm)

  W = W_matrix(ps,pm)(γαβ_s)
  u_m = inv_jacobian(pm)(γαβ_m) ⋅ fV(pm)(γαβ_m)
  u_s = W ⋅ u_m
  tangent_s = W ⋅ tangent_m
  u_s ⋅ tangent_s
end

dc_s = ∫(  u_dot_t_p1   )dΩ_edge
dc_s.dict
sum(dc_s) ≈ sum(dc_m)
sum(dc_s) != dof_1

;
################################################################################
########### Edge [9 11] #######################
################################################################################

## Panel 4 is the master, panel 2 is the slave
pm = 4
ps = 2

## In panel 4, edge [9 11] -> dof 10
## In panel 2, edge [9 11] -> dof 6
dof_4 = cell_dofs[4][10]
dof_2 = cell_dofs[2][6]


############# PANEL 4
## In panel 4, edge [9 11] -> (γ,α,β) = (0,π/4,β) where β = [-π/4,π/4], tangent = (0,0,1)
dof_4
tangent_m = VectorValue(0.0,0.0,1.0)
function u_dot_t_p4(β)
  γαβ = VectorValue(0.0,π/4,β[1])
  u = inv_jacobian(pm)(γαβ) ⋅ fV(pm)(γαβ)
  u ⋅ tangent_m
end
dc_m = ∫(  u_dot_t_p4   )dΩ_edge
dc_m.dict
sum(dc_m) ≈ dof_4 ### Why is this false?

############# PANEL 2
## In panel 2, edge [9 11] -> (γ,α,β) = (0,α,π/4) where α = [-π/4,π/4], tangent = (0,1,0)
dof_2
function u_dot_t_p2(α)
  γαβ_s = VectorValue(0.0,α[1],π/4)
  γαβ_m = slave2master(γαβ_s,ps,pm)

  W = W_matrix(ps,pm)(γαβ_s)
  u_m = inv_jacobian(pm)(γαβ_m) ⋅ fV(pm)(γαβ_m)
  u_s = W ⋅ u_m
  tangent_s = W ⋅ tangent_m
  u_s ⋅ tangent_s
end

dc_s = ∫(  u_dot_t_p2   )dΩ_edge
dc_s.dict
sum(dc_s) ≈ sum(dc_m)
sum(dc_s) != dof_2



################################################################################
########### Edge [7 11] #######################
################################################################################

########### Panel 3 is the master, panel 2 is the slave
pm = 3
ps = 2

## In panel 3, edge [7 11] -> dof 6
## In panel 2, edge [7 11] -> dof 10
dof_3 = cell_dofs[3][6]
dof_2 = cell_dofs[2][10]

############# PANEL 3
## In panel 3, edge [7 11] -> (γ,α,β) = (0,α,π/4) where α = [-π/4,π/4], tangent = (0,1,0)
dof_3
tangent_m = VectorValue(0.0,1.0,0.0)
function u_dot_t_p3(α)
  γαβ = VectorValue(0.0,α[1],π/4)
  u = inv_jacobian(pm)(γαβ) ⋅ fV(pm)(γαβ)
  u ⋅ tangent_m
end
dc_m = ∫(  u_dot_t_p3   )dΩ_edge
dc_m.dict
sum(dc_m) ≈ dof_3

############# PANEL 2
## In panel 2, edge [7 11] -> (γ,α,β) = (0,π/4,β) where β = [-π/4,π/4], tangent = (0,0,1)
dof_2
function u_dot_t_p2(β)
  γαβ_s = VectorValue(0.0,π/4,β[1])
  γαβ_m = slave2master(γαβ_s,ps,pm)

  W = W_matrix(ps,pm)(γαβ_s)
  u_m = inv_jacobian(pm)(γαβ_m) ⋅ fV(pm)(γαβ_m)
  u_s = W ⋅ u_m
  tangent_s = W ⋅ tangent_m
  u_s ⋅ tangent_s
end

dc_s = ∫(  u_dot_t_p2   )dΩ_edge
dc_s.dict
sum(dc_s) ≈ sum(dc_m)
