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
using BenchmarkTools
include("cube_cell_maps.jl")

struct ManifoldGrid{Dc,Dp,Dp_topo,Tp,O,Tn,A} <: Grid{Dc,Dp} where Tn
  cube_grid::UnstructuredGrid{Dc,Dp_topo,Tp,O,Tn}
  phys_grid::UnstructuredGrid{Dc,Dp,Tp,O,Tn} ## grid where integration occurs
  phys_cell_map::A ## map to integration space
  # manifold_cell_map::A ### map to manifold - not the integration space
end

function ManifoldGrid(cube_model::DiscreteModel)
  cube_grid = get_grid(cube_model)
  panel_ids = get_panel_ids(cube_model)
  ManifoldGrid(cube_grid,panel_ids)
end

function ManifoldGrid(cube_grid::UnstructuredGrid{Dc,Dp_topo,Tp,O,Tn},panel_ids::AbstractArray{<:Int64}) where {Dc,Dp_topo,Tp,O,Tn}

  cell_coords = get_cell_coordinates(cube_grid)
  cell_node_ids = get_cell_node_ids(cube_grid)
  cmaps = get_cell_map(cube_grid)

  ## map the cube coordinates -> physical integration space
  coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
  coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)
  coords_panel1_3D = lazy_map(BumpMap(), coords_panel1_2D)
  coords_panelp = lazy_map(InvPanelMap(), coords_panel1_3D, panel_ids)

  ## evaluate the nodes
  T = eltype(eltype(coords_panelp))
  phys_nodes = similar(coords_panelp, T, num_nodes(cube_grid))
  get_nodes_from_coords!(phys_nodes,cell_node_ids,coords_panelp)

  ## define grid of integratin space
  phys_grid = Gridap.Geometry.UnstructuredGrid(phys_nodes,cell_node_ids,
      get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid))
  Dp = num_point_dims(phys_grid)

  ## cell map to integration
  phys_cell_map = lazy_map(CellPanelMaps(Rp1,R1p,Bump), panel_ids, cmaps)
  A = typeof(phys_cell_map)


  ManifoldGrid{Dc,Dp,Dp_topo,Tp,O,Tn,A}(cube_grid,phys_grid,phys_cell_map)

end


"""
grid API
"""
function Gridap.Geometry.get_node_coordinates(grid::ManifoldGrid)
  Gridap.Geometry.get_node_coordinates(grid.phys_grid)
end

function Gridap.Geometry.get_cell_node_ids(grid::ManifoldGrid)
  Gridap.Geometry.get_cell_node_ids(grid.phys_grid)
end

function Gridap.Geometry.get_reffes(grid::ManifoldGrid)
  Gridap.Geometry.get_reffes(grid.phys_grid)
end

function Gridap.Geometry.get_cell_type(grid::ManifoldGrid)
  Gridap.Geometry.get_cell_type(grid.phys_grid)
end

function Gridap.Geometry.get_cell_map(grid::ManifoldGrid)
  grid.phys_cell_map
end

function get_manifold_cell_map(grid::ManifoldGrid)
  grid.manifold_cell_map
end

Gridap.Geometry.OrientationStyle(grid::ManifoldGrid) = Gridap.Geometry.NonOriented()


####
PanelMap() = PanelRotationMap(rp1)
InvPanelMap() = PanelRotationMap(r1p)
BumpMap() = Panel1BumpMap(A_bump,B_bump,b_bump)

phys_grid = ManifoldGrid(cube_model_3D)
phys_cmap = get_cell_map(phys_grid)

test_cell_maps(phys_cmap,ref_cell_coords,cell_coords)


ref_model = Gridap.Adaptivity.refine(cube_model_3D)

ref_grid = ManifoldGrid(ref_model)
test_cell_maps(get_cell_map(ref_grid),get_cell_ref_coordinates(ref_model),get_cell_coordinates(ref_model))

# writevtk(ref_grid,dir*"/ref_grid",append=false)
