using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Test
using LinearAlgebra
using FillArrays
include("maps/panel_rotations.jl")
include("maps/bump_panel1.jl")
include("maps/maps.jl")
include("maps/panel_ids_from_refinement_v2.jl")
include("cube_topo/cube_surface_1_cell_per_panel.jl")


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

cube_model_3D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel()...)

_A , _B, _b = bump_matrics()
A_bump = TensorValue(_A)
B_bump = TensorValue(_B)
b_bump = VectorValue(_b)

BumpMap() = Panel1BumpMap(A_bump,B_bump,b_bump)

#####

model = cube_model_3D
ref_model = Gridap.Adaptivity.refine(model)
ref_ref_model = Gridap.Adaptivity.refine(ref_model)
ref_ref_ref_model = Gridap.Adaptivity.refine(ref_ref_model)

# coordinates of panel 1
panel_ids = get_panel_ids(ref_ref_ref_model)
cell_coords = get_cell_coordinates(ref_ref_ref_model)
p1 = findall(x->x==1,panel_ids)
cell_coords_panel1 = cell_coords[p1]
cell_coords_panel1_2D = lazy_map(BumpMap(), cell_coords_panel1)

### a big hack is to refine a cartesian panel 1 at the same time as the cube model
### then the reordering is not needed ....
### seems very dirty to me, but it works
panel1 = UnstructuredDiscreteModel(CartesianGrid((-1,1,-1,1), (1,1) ))
ref_panel1 = Gridap.Adaptivity.refine(panel1)
ref_ref_panel1 = Gridap.Adaptivity.refine(ref_panel1)
ref_ref_ref_panel1 = Gridap.Adaptivity.refine(ref_ref_panel1)

_cell_coords_panel1 = get_cell_coordinates(ref_ref_ref_panel1)

println(cell_coords_panel1_2D .== _cell_coords_panel1)


################################################################################

# coordinates of panel 1
panel_ids = get_panel_ids(ref_ref_ref_model)
cell_coords = get_cell_coordinates(ref_ref_ref_model)
p1 = findall(x->x==1,panel_ids)
cell_coords_panel1 = cell_coords[p1]
cell_coords_panel1_2D = lazy_map(BumpMap(), cell_coords_panel1)

## separate grid for panel 1
nC_fine = Tuple(fill(Int(sqrt(length(p1))),2))
fine_panel1 = UnstructuredGrid(CartesianGrid((-1,1,-1,1), nC_fine ))
_cell_coords_panel1 = get_cell_coordinates(fine_panel1)

shuffled_panel1_cell_coords = copy(_cell_coords_panel1)

#### START
model = ref_ref_ref_model
parent = model.parent

Nnew =  Int( num_cells(model)/6 )
Nold = Int(  num_cells(parent)/6 )
v = Tuple( fill( Int(sqrt(Nnew) ), 2) )
w = Tuple( fill(Int(sqrt(Nold)  ), 2) )
z = (2,2)

f2c_cell_map, fcell_to_child_id = Gridap.Adaptivity._create_cartesian_f2c_maps((4,4), (2,2))
c2f_cells = Adaptivity.get_o2n_faces_map(f2c_cell_map)
reindex = c2f_cells.data

target_ids  = reduce(vcat, [reindex[1:8], reindex[17:24],
                            reindex[9:16], reindex[25:32],
                            reindex[33:40], reindex[49:56],
                            reindex[41:48], reindex[57:64]])
println(cell_coords_panel1_2D .== _cell_coords_panel1[target_ids])

_f2c_cell_map,  = Gridap.Adaptivity._create_cartesian_f2c_maps((2,2), (2,2))

o_f2c_cell_map,  = Gridap.Adaptivity._create_cartesian_f2c_maps((1,1), (2,2))


ids = findall(x->x==1,_f2c_cell_map)

new = []

#### start recusion
glue = get_adaptivity_glue(model)

test = (glue.n2o_faces_map[end][p1])


#####


_glue = get_adaptivity_glue(model.parent)

_test = (_glue.n2o_faces_map[end][p1])


_c2f_cells = Adaptivity.get_o2n_faces_map(_f2c_cell_map)
reindex = c2f_cells.data



###### RECURSE 1
model = ref_ref_model
glue = get_adaptivity_glue(model)
_test = (glue.n2o_faces_map[end][p1])


###### RECURSE 2
model = ref_model
glue = get_adaptivity_glue(model)
o_test = (glue.n2o_faces_map[end][p1])

error

### finish -- do nothing

#

# ################################################################################
# ##### RECURSE 1
# model = ref_ref_model

# parent = model.parent

# Nnew =  Int( num_cells(model)/6 )
# Nold = Int(  num_cells(parent)/6 )
# v = Tuple( fill( Int(sqrt(Nnew) ), 2) )
# w = Tuple( fill(Int(sqrt(Nold)  ), 2) )
# z = (4,4)


# f2c_cell_map, fcell_to_child_id = Gridap.Adaptivity._create_cartesian_f2c_maps((2,2),(4,4))
# c2f_cells = Adaptivity.get_o2n_faces_map(f2c_cell_map)
# _reindex = c2f_cells.data



# println(cell_coords_panel1_2D .== shuffled_panel1_cell_coords[_reindex[reindex]])

# ###### RECURSE 2
# model = ref_model
# parent = model.parent

# Nnew =  Int( num_cells(model)/6 )
# Nold = Int(  num_cells(parent)/6 )
# v = Tuple( fill( Int(sqrt(Nnew) ), 2) )
# w = Tuple( fill(Int(sqrt(Nold)  ), 2) )
# z = (8,8)

# f2c_cell_map, fcell_to_child_id = Gridap.Adaptivity._create_cartesian_f2c_maps((1,1),(8,8))
# c2f_cells = Adaptivity.get_o2n_faces_map(f2c_cell_map)
# final_reindex = c2f_cells.data

# println(cell_coords_panel1_2D .== shuffled_panel1_cell_coords[reindex[_reindex[final_reindex]]])


# ####
