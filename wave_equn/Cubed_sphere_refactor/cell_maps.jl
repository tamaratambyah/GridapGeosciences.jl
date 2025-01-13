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


ncells_per_panel = 4
npanels = 6
ref_panel = CartesianDiscreteModel((-1,1,-1,1),(ncells_per_panel,ncells_per_panel))

cell_coordinates = get_cell_coordinates(ref_panel)
cell_node_ids = get_cell_node_ids(ref_panel)
node_coordinates = get_node_coordinates(ref_panel)

panel1 = get_panel(1,node_coordinates,cell_node_ids)

cmaps = get_cell_map(panel1)./1
evaluate(cmaps[1],Point(1,1))
