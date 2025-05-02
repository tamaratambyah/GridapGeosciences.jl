using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../src/initialise.jl")

nodes = [
  Point(1.0, -1.0, -1.0)  # node 1
  Point(1.0, 1.0, -1.0)   # node 2
  Point(1.0, -1.0, 1.0)   # node 3
  Point(1.0, 1.0, 1.0)    # node 4
]


Rp1, R1p = panel_rotations()

for i in 1:6
  println("Panel $i")
  p_nodes = R1p[i] .⋅ nodes
  println(p_nodes)

  _nodes = Rp1[i] .⋅ p_nodes
  @test _nodes == nodes
  println(_nodes)
end

@test Rp1 == inv.(R1p)
@test R1p == inv.(Rp1)


coarse_model = ManifoldDiscreteModel(cube_model_3D,cube)
panel_ids = get_panel_ids(coarse_model)

parametric_grid = get_parametric_grid(get_grid((coarse_model)))
ambient_grid = get_ambient_grid(get_grid(coarse_model))

ambient_nodes = get_node_coordinates(ambient_grid)
cell_coords = get_cell_coordinates(ambient_grid)

lazy_map(PanelRotationField(rp1),ambient_nodes,panel_ids)
lazy_map(PanelRotationMap(rp1),ambient_nodes,panel_ids)

lazy_map(PanelRotationMap(rp1), cell_coords, panel_ids)
lazy_map(PanelRotationField(rp1), cell_coords, panel_ids)




coarse_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
model = (Gridap.Adaptivity.refine(coarse_model))
panel_ids = get_panel_ids(model)
ambient_grid = get_ambient_grid(get_grid(model))
cell_coords = get_cell_coordinates(ambient_grid)
cell_node_ids = get_cell_node_ids(ambient_grid)

_cell_coords = lazy_map(Sigma(), cell_coords)

plot_coords(_cell_coords,panel_ids,cell_node_ids)
