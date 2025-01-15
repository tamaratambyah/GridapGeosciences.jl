using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Test
using LinearAlgebra
using FillArrays
include("individual_panels.jl")


ncells_per_panel = 2
npanels = 6
ref_panel = CartesianDiscreteModel((-1,1,-1,1),(ncells_per_panel,ncells_per_panel))
ref_cube = CartesianDiscreteModel((-1,1,-1,1,-1,1),(ncells_per_panel,ncells_per_panel,ncells_per_panel))


cell_coordinates = get_cell_coordinates(ref_panel)
cell_node_ids = get_cell_node_ids(ref_panel)
node_coordinates = get_node_coordinates(ref_panel)

# append first 2 panels
panel1 = get_panel(1,node_coordinates,cell_node_ids)
panel2 = get_panel(2,node_coordinates,cell_node_ids)
append_panels = lazy_append(panel1,panel2)

# append panels 3--6
for panel = collect(3:npanels)
  cube_panel = get_panel(panel,node_coordinates,cell_node_ids)
  append_panels = lazy_append(append_panels,cube_panel)
end


cubed_sphere_grid = UnstructuredGrid(append_panels)
writevtk(cubed_sphere_grid,datadir("CubedSphere")*"/cube_sphere_grid",append=false)
writevtk(ref_cube,datadir("CubedSphere")*"/ref_cube",append=false)

cube_nodes = sort(get_node_coordinates(cubed_sphere_grid))
cube_cell_node_ids = get_cell_node_ids(cubed_sphere_grid)

### Append does not remove the duplicates, I need to do this manually ....
