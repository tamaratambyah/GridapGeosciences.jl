using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson
using Gridap.Geometry

include("../convergence_tools.jl")

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


dir = datadir("InterpolationConvergence_glued")
!isdir(dir) && mkdir(dir)

models = get_octree_refined_models(ranks,1,true)
panel_model = models[end].models.item

panel_ids = get_panel_ids(panel_model)

grid = get_grid(panel_model)
topo = get_grid_topology(panel_model)
cmap = get_cell_map(grid)

using Gridap.Geometry
Gridap.Geometry.get_faces(topo,2,0)

coords = QUAD.vertex_coords

masters = [6, 5, 6, 3, 6, 4, 5, 6]
m_vertex = [1,3,2,3,4,4,4,3]
slaves = [[1,5], [1,3], [1,2], [1,2], [2,4], [2,3], [3,4], [4,5]]
s_vertex = [ [1,1], [2,1], [3,1], [4,2], [3,2], [4,4], [2,3], [1,2]  ]

i = 1

pm = masters[i]
xm = evaluate(cmap[pm],coords[m_vertex[i]])
Jm = forward_jacobian(pm)(xm)

ps1 = slaves[i][1]
xs1 = evaluate(cmap[ps1],coords[s_vertex[i][1]])
Jinvs1 = forward_pinv_jacobian(ps1)(xs1)
coeffs1 = Matrix(Jinvs1⋅Jm)

ps2 = slaves[i][2]
xs2 = evaluate(cmap[ps2],coords[s_vertex[i][2]])
Jinvs2 = forward_pinv_jacobian(ps2)(xs2)
coeffs2 = Matrix(Jinvs2⋅Jm)
