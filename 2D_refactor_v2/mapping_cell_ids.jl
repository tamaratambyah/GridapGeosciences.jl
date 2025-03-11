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
_model = ref_ref_model
panel_ids = get_panel_ids(_model)
cell_coords = get_cell_coordinates(_model)
p1 = findall(x->x==1,panel_ids)
cell_coords_panel1 = cell_coords[p1]
cell_coords_panel1_2D = lazy_map(BumpMap(), cell_coords_panel1)

### a big hack is to refine a cartesian panel 1 at the same time as the cube model
### then the reordering is not needed ....
### seems very dirty to me, but it works
# panel1 = UnstructuredDiscreteModel(CartesianGrid((-1,1,-1,1), (1,1) ))
# ref_panel1 = Gridap.Adaptivity.refine(panel1)
# ref_ref_panel1 = Gridap.Adaptivity.refine(ref_panel1)
# ref_ref_ref_panel1 = Gridap.Adaptivity.refine(ref_ref_panel1)
# _cell_coords_panel1 = get_cell_coordinates(ref_ref_ref_panel1)
# println(cell_coords_panel1_2D .== _cell_coords_panel1)


################################################################################
num_cells_p1 = length(p1)

## create cartesian grid for panel 1
nC_p1 = Tuple(fill(Int(sqrt(num_cells_p1)),2))
panel1 = UnstructuredGrid(CartesianGrid((-1,1,-1,1), nC_p1 ))
_cell_coords_panel1 = get_cell_coordinates(panel1)

## now do the reordering
n = Int( num_cells_p1/4) # 1:4 partition rule
m = Int(n/4)

nc_fine = Tuple(fill(Int(sqrt(n)),2))
nc_coarse = Tuple(fill(Int(sqrt(m)),2))

f2c_cell_map, fcell_to_child_id = Gridap.Adaptivity._create_cartesian_f2c_maps(nc_fine, (2,2))
c2f_cell_map = Adaptivity.get_o2n_faces_map(f2c_cell_map)

_f2c_cell_map,  = Gridap.Adaptivity._create_cartesian_f2c_maps(nc_coarse, (2,2))

p = sortperm(_f2c_cell_map)

reindex = c2f_cell_map[p].data

println(cell_coords_panel1_2D .== _cell_coords_panel1[reindex])
