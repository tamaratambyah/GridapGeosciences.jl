using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using Gridap.FESpaces, Gridap.ReferenceFEs
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
# include("../../Geophysical/CurlConformingFESpacesFixes.jl")

## pullback 3D vector to 3D chart
inv_jacobian(p) = x -> inv(forward_jacobian(p)(x))
contra_v_3D(vecX::Function,p) = αβ -> inv_jacobian(p)(αβ) ⋅ vecX(p)(αβ)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)

transpose_jacobian(p) = x -> transpose(forward_jacobian(p)(x))
inv_tranpose_jacobian(p) = x -> inv(transpose_jacobian(p)(x))
covar_v_3D(vecX::Function,p) = αβ -> transpose_jacobian(p)(αβ) ⋅ vecX(p)(αβ)
covar_v_3D(vecX::Function) = p -> covar_v_3D(vecX,p)

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

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
### pull with JT -> to a 1-form i.e. covariant vector (index down)
vec_cov_cf = ParametricCellField(covar_v_3D(fV),Ω,panel_ids)

### lowest order nedelec
order = 0
et = Float64
reffe =  ReferenceFE(nedelec,et,order)
R = TestFESpace(panel_model,reffe;conformity=:Hcurl,dirichlet_tags=tags)
H = TrialFESpace(R,vec_cov_cf)

dof_basis = get_fe_dof_basis(H)
vec_cov_h = dof_basis(vec_cov_cf)
cell_dofs = collect(vec_cov_h)

#### Gamma = 0 edges:
# edge [3 7]:
cell_dofs[1][11] # Panel 1 dof 11
cell_dofs[3][9] # panl 3 dof 9

# edge [3 13]:
cell_dofs[3][5] # Panel 3 dof 5
cell_dofs[5][7] # panel 5 dof 7

# edge [1 3]:
cell_dofs[1][5]
cell_dofs[5][9]

# edge [1 5]:
cell_dofs[1][9]
cell_dofs[6][5]

# edge [5 9]:
cell_dofs[2][9]
cell_dofs[6][11]

# edge [15 9]:
cell_dofs[4][5]
cell_dofs[6][7]

# edge [9 11]:
cell_dofs[2][7]
cell_dofs[4][11]

# edge [5 7]:
cell_dofs[1][7]
cell_dofs[2][5]

# edge [7 11]:
cell_dofs[2][11]
cell_dofs[3][7]

# edge [13 11]:
cell_dofs[3][11]
cell_dofs[4][7]

# edge [15 13]:
cell_dofs[4][9]
cell_dofs[5][11]

# edge [1 15]:
cell_dofs[5][5]
cell_dofs[6][9]



#### Gamma = 1 edges:
# edge [4 8]:
cell_dofs[1][12]
cell_dofs[3][10]

# edge [4 14]:
cell_dofs[3][6] # Panel 3 dof 5
cell_dofs[5][8] # panel 5 dof 7

# edge [2 4]:
cell_dofs[1][6]
cell_dofs[5][10]

# edge [2 6]:
cell_dofs[1][10]
cell_dofs[6][6]

# edge [6 10]:
cell_dofs[2][10]
cell_dofs[6][12]

# edge [16 10]:
cell_dofs[4][6]
cell_dofs[6][8]

# edge [10 12]:
cell_dofs[2][8]
cell_dofs[4][12]

# edge [6 8]:
cell_dofs[1][8]
cell_dofs[2][6]

# edge [8 12]:
cell_dofs[2][12]
cell_dofs[3][8]

# edge [14 12]:
cell_dofs[3][12]
cell_dofs[4][8]

# edge [16 14]:
cell_dofs[4][10]
cell_dofs[5][12]

# edge [2 16]:
cell_dofs[5][6]
cell_dofs[6][10]


#### Vertical edges, variable gamma:
# edge [1 2]:
cell_dofs[1][1]
cell_dofs[5][1]
cell_dofs[6][1]

# edge [3 4]:
cell_dofs[1][2]
cell_dofs[3][1]
cell_dofs[5][3]

# edge [5 6]:
cell_dofs[1][3]
cell_dofs[2][1]
cell_dofs[6][2]

# edge [7 8]:
cell_dofs[2][2]
cell_dofs[3][3]
cell_dofs[1][4]

# edge [9 10]:
cell_dofs[2][3]
cell_dofs[4][2]
cell_dofs[6][4]

# edge [11 12]:
cell_dofs[2][4]
cell_dofs[3][4]
cell_dofs[4][4]

# edge [13 14]:
cell_dofs[3][2]
cell_dofs[4][3]
cell_dofs[5][4]

# edge [15 16]:
cell_dofs[4][1]
cell_dofs[5][2]
cell_dofs[6][3]



# topo = get_grid_topology(panel_model)
# Dc = num_cell_dims(topo)
# edge2cell  = Gridap.Geometry.get_faces(topo, 1, Dc)
# cells2edge  = Gridap.Geometry.get_faces(topo, Dc, 1)
# edge2node = get_faces(topo,1,0)

# tangents = get_edge_tangent(HEX)


# cell = 1

# edge_ids = cells2edge[cell]

# edges = edge2cell[edge_ids]

# edge = edges[11]
# pm = edge[end]
# ps = cell

# 11;
# ################################################################################
# ########### Edge [3 7] #######################
# ################################################################################

# ########### Panel 3 is the master, panel 1 is the slave
# pm = 3
# ps = 1

# ## In panel 3,edge [3 7] -> dof 9
# ## In panel 1, edge [3 7] -> dof 10
# dof_3 = cell_dofs[3][9]
# dof_1 = cell_dofs[1][10]

# ## In panel 3, Node 3 -> (γ,α,β) = (0,-π/4,-π/4)
# ## In panel 1, Node 3 -> (γ,α,β) = (0,π/4,-π/4)
# vm = Point(0,-π/4,-π/4)
# vs = Point(0,π/4,-π/4)

# lag_reffe = Gridap.ReferenceFEs.LagrangianRefFE(Float64,SEGMENT,1)
# shapefuns = get_shapefuns(lag_reffe)
# edge_map = linear_combination([vm, Point(0,-π/4,π/4)], shapefuns)
# jac_edge_map = gradient(edge_map)
# measure = meas(jac_edge_map)

# pt = [vm, Point(0,-π/4,π/4)]
# measure(Point(1))

# function slave2master(γαβ_s,p_s,p_m)
#   ## 1. push slave to ambient space
#   ## 2. pull to master using inverse of master map
#   xyz = ForwardMap(p_s)(γαβ_s)
#   inv_m = inverse_map(ForwardMap(p_m))
#   γαβ_m = evaluate(inv_m,xyz)
#   return γαβ_m
# end

# #### Node 3: given the slave parametric coord, compute the master parametric coord
# γαβ_s = Point(0,π/4,-π/4)
# γαβ_m = slave2master(γαβ_s,ps,pm)
# γαβ_m ≈ vm


# #### Tangent to the edge [3, 7]
# tangent_m = VectorValue(0.0,0.0,1.0)

# W = inv_jacobian(ps)(γαβ_s)⋅forward_jacobian(pm)(γαβ_m)
# coeffs = Matrix(W)
# tangent_s = W ⋅ tangent_m
# tangent_s ≈ VectorValue(0,0,1)


# #### Cofficient matrix in terms of slave points
# function W_matrix(p_s,p_m)
#   function _W_matrix(γαβ_s)
#     W = inv_jacobian(p_s)(γαβ_s)⋅forward_jacobian(p_m)( slave2master(γαβ_s,p_s,p_m) )
#     W
#   end
# end

# W_s = W_matrix(ps,pm)(γαβ_s)
# coeffs_s = Matrix(W_s)
# coeffs_s ≈ coeffs


# ########### integral over the edge:
# ##### an the edge [3,7] is [0,-π/4,-π/4] -> [0,π/4,π/4]
# ##### we could restrict the 3D model to the edge, but as a hack, let's create a
# ##### 1D model that is the line [-π/4,π/4]
# model_edge = UnstructuredDiscreteModel(CartesianDiscreteModel((-π/4,π/4),(1)))
# Ω_edge = Triangulation(model_edge)
# dΩ_edge = Measure(Ω_edge,0) ## zero order elements -> zero order quadrature for edge values


# ############# PANEL 3
# ## In panel 3,edge [3 7] -> (γ,α,β) = (0,-π/4,β) where β = [-π/4,π/4], tangent = (0,0,1)
# tangent_m = VectorValue(0,0,1)
# function u_dot_t_p3(β)
#   γαβ = VectorValue(0.0,-π/4,β[1])
#   u = inv_jacobian(pm)(γαβ) ⋅ fV(pm)(γαβ)
#   u ⋅ tangent_m
# end
# dc_m = ∫(  u_dot_t_p3   )dΩ_edge
# dc_m.dict
# sum(dc_m) ≈ dof_3


# ############# PANEL 1
# ## In panel 1, edge [3 7] -> (γ,α,β) = (0,π/4,β) where β = [-π/4,π/4], tangent = (0,0,1)
# dof_1
# function u_dot_t_p1(β)
#   γαβ_s = VectorValue(0.0,π/4,β[1])
#   γαβ_m = slave2master(γαβ_s,ps,pm)

#   W = W_matrix(ps,pm)(γαβ_s)
#   u_m = inv_jacobian(pm)(γαβ_m) ⋅ fV(pm)(γαβ_m)
#   u_s = W ⋅ u_m
#   tangent_s = W ⋅ tangent_m
#   u_s ⋅ tangent_s
# end

# dc_s = ∫(  u_dot_t_p1   )dΩ_edge
# dc_s.dict
# sum(dc_s) ≈ sum(dc_m)
# sum(dc_s) != dof_1

# ;
# ################################################################################
# ########### Edge [9 11] #######################
# ################################################################################

# ## Panel 4 is the master, panel 2 is the slave
# pm = 4
# ps = 2

# ## In panel 4, edge [9 11] -> dof 10
# ## In panel 2, edge [9 11] -> dof 6
# dof_4 = cell_dofs[4][10]
# dof_2 = cell_dofs[2][6]


# ############# PANEL 4
# ## In panel 4, edge [9 11] -> (γ,α,β) = (0,π/4,β) where β = [-π/4,π/4], tangent = (0,0,1)
# dof_4
# tangent_m = VectorValue(0.0,0.0,1.0)
# function u_dot_t_p4(β)
#   γαβ = VectorValue(0.0,π/4,β[1])
#   u = inv_jacobian(pm)(γαβ) ⋅ fV(pm)(γαβ)
#   u ⋅ tangent_m
# end
# dc_m = ∫(  u_dot_t_p4   )dΩ_edge
# dc_m.dict
# sum(dc_m) ≈ dof_4 ### Why is this false?

# ############# PANEL 2
# ## In panel 2, edge [9 11] -> (γ,α,β) = (0,α,π/4) where α = [-π/4,π/4], tangent = (0,1,0)
# dof_2
# function u_dot_t_p2(α)
#   γαβ_s = VectorValue(0.0,α[1],π/4)
#   γαβ_m = slave2master(γαβ_s,ps,pm)

#   W = W_matrix(ps,pm)(γαβ_s)
#   u_m = inv_jacobian(pm)(γαβ_m) ⋅ fV(pm)(γαβ_m)
#   u_s = W ⋅ u_m
#   tangent_s = W ⋅ tangent_m
#   u_s ⋅ tangent_s
# end

# dc_s = ∫(  u_dot_t_p2   )dΩ_edge
# dc_s.dict
# sum(dc_s) ≈ sum(dc_m)
# sum(dc_s) != dof_2



# ################################################################################
# ########### Edge [7 11] #######################
# ################################################################################

# ########### Panel 3 is the master, panel 2 is the slave
# pm = 3
# ps = 2

# ## In panel 3, edge [7 11] -> dof 6
# ## In panel 2, edge [7 11] -> dof 10
# dof_3 = cell_dofs[3][6]
# dof_2 = cell_dofs[2][10]

# ############# PANEL 3
# ## In panel 3, edge [7 11] -> (γ,α,β) = (0,α,π/4) where α = [-π/4,π/4], tangent = (0,1,0)
# dof_3
# tangent_m = VectorValue(0.0,1.0,0.0)
# function u_dot_t_p3(α)
#   γαβ = VectorValue(0.0,α[1],π/4)
#   u = inv_jacobian(pm)(γαβ) ⋅ fV(pm)(γαβ)
#   u ⋅ tangent_m
# end
# dc_m = ∫(  u_dot_t_p3   )dΩ_edge
# dc_m.dict
# sum(dc_m) ≈ dof_3

# ############# PANEL 2
# ## In panel 2, edge [7 11] -> (γ,α,β) = (0,π/4,β) where β = [-π/4,π/4], tangent = (0,0,1)
# dof_2
# function u_dot_t_p2(β)
#   γαβ_s = VectorValue(0.0,π/4,β[1])
#   γαβ_m = slave2master(γαβ_s,ps,pm)

#   W = W_matrix(ps,pm)(γαβ_s)
#   u_m = inv_jacobian(pm)(γαβ_m) ⋅ fV(pm)(γαβ_m)
#   u_s = W ⋅ u_m
#   tangent_s = W ⋅ tangent_m
#   u_s ⋅ tangent_s
# end

# dc_s = ∫(  u_dot_t_p2   )dΩ_edge
# dc_s.dict
# sum(dc_s) ≈ sum(dc_m)
