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
include("initialise.jl")
# include("cube_cell_maps.jl")
# include("sphere_nodes_via_panel1.jl")
# include("cube_via_panel1.jl")

struct ManifoldGrid{Dc,Dp,Dp_topo,Dp_phys,A<:Grid{Dc,Dp_topo},B<:Grid{Dc,Dp_phys},C<:Grid{Dc,Dp},E} <: Grid{Dc,Dp}
  cube_grid::A
  phys_grid::B ## physical domain (for FEM)
  ambient_grid::C ## ambient space (manifold)
  phys_cell_map::E ## map to integration space
  # manifold_cell_map::A ### map to manifold - not the integration space
end

function ManifoldGrid(cube_grid<:Grid{Dc,Dp_topo},phys_grid<:Grid{Dc,Dp_phys},
  ambient_grid<:Grid{Dc,Dp},phys_cell_map) where {Dc,Dp,Dp_topo,Dp_phys}

  A = typeof(cube_grid)
  B = typeof(phys_grid)
  C = typeof(ambient_grid)
  E = typeof(phys_cell_map)
  ManifoldGrid{Dc,Dp,Dp_topo,Dp_phys,A,B,C,E}(cube_grid,phys_grid,ambient_grid,phys_cell_map)
end




function CSGrid(cube_model::DiscreteModel,nodes_function::Function,CellPanelMap::Map)
  cube_grid = get_grid(cube_model)
  panel_ids = get_panel_ids(cube_model)
  phys_coords,  manifold_coords,  = nodes_function(cube_model)
  ManifoldGrid(cube_grid,panel_ids,phys_coords,manifold_coords,CellPanelMap)
end

function ManifoldGrid(cube_grid::UnstructuredGrid{Dc,Dp_topo,Tp,O,Tn},
  panel_ids::AbstractArray{<:Int64},phys_coords,manifold_coords,
  CellPanelMap::Map) where {Dc,Dp_topo,Tp,O,Tn}

  cell_coords = get_cell_coordinates(cube_grid)
  cell_node_ids = get_cell_node_ids(cube_grid)
  cmaps = get_cell_map(cube_grid)

  ## evaluate the physical nodes
  phys_nodes = get_nodes_from_coords(cube_grid,phys_coords)

  ## define grid of physical space (for FEM)
  phys_grid = Gridap.Geometry.UnstructuredGrid(phys_nodes,cell_node_ids,
      get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid))
  Dp = num_point_dims(phys_grid)

  ## cell map to integration
  phys_cell_map = lazy_map(CellPanelMap, panel_ids, cmaps)
  A = typeof(phys_cell_map)


  ## evaluate the manifold nodes
  manifold_nodes = get_nodes_from_coords(cube_grid,manifold_coords)

  ## define grid of ambient space (manifold)
  manifold_grid = Gridap.Geometry.UnstructuredGrid(manifold_nodes,cell_node_ids,
      get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid))



  ManifoldGrid{Dc,Dp,Dp_topo,Tp,O,Tn,A}(cube_grid,phys_grid,phys_cell_map,manifold_grid)

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

function get_manifold_grid(grid::ManifoldGrid)
  grid.manifold_grid
end

function get_cube_grid(grid::ManifoldGrid)
  grid.cube_grid
end

Gridap.Geometry.OrientationStyle(grid::ManifoldGrid) = Gridap.Geometry.NonOriented()


####
_model = ref_ref_ref_model

### make the cube
phys_grid = CSGrid(_model,cube_nodes,CubeCellPanelMaps(Rp1,R1p,Bump))
writevtk(get_manifold_grid(phys_grid),dir*"/C_grid",append=false)

phys_cmap = get_cell_map(phys_grid)
test_cell_maps(phys_cmap, get_cell_ref_coordinates(get_grid(_model)), get_cell_coordinates(_model))

### make the cubed sphere
phys_grid = CSGrid(_model,cubed_sphere_nodes,CellPanelMaps(Rp1,R1p,Bump))
writevtk(get_manifold_grid(phys_grid),dir*"/CS_grid",append=false)

phys_cmap = get_cell_map(phys_grid)
test_cell_maps(phys_cmap, get_cell_ref_coordinates(get_grid(_model)), get_cell_coordinates(_model))
