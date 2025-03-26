"""
    ManifoldGrid

Stores the information required to solve on the surface

topo_grid               : grid of the topology (i.e. cube)
parametric_grid         : grid of the parametric domain (for FEM) i.e. panel 1
                          - note, this grid represents panel 1, and only panel 1
                          nodes are non-zero.
ambient_grid            : grid of the ambient space (manifold)
parametric_cell_map     : map to the parametric space (used for integration)
ambient_cell_map        : map to ambient space (manifold - not the integration space)
parametric_cell_coords  : the cell coords in parametric space
                        - this must be stored separatately as we cannot have dG
                         type nodes
panel_ids               : panel ids corresponding to the panels of the cube

The grid interface is overloaded to return the parametric_grid
- [`get_node_coordinates(grid::ManifoldGrid)]
- [`get_cell_node_ids(grid::ManifoldGrid)`]
- [`get_reffes(grid::ManifoldGrid)`]
- [`get_cell_type(grid::ManifoldGrid)`]
- [`get_cell_ref_coordinates(grid::ManifoldGrid)`]
- [`get_cell_coordinates(grid::ManifoldGrid)`]
- ['OrientationStyle(grid::ManifoldGrid)``]

The following methods have been added to return the ambient_grid
- [`get_ambient_node_coordinates(grid::ManifoldGrid)`]
- [`get_ambient_cell_coordinates(grid::ManifoldGrid)`]
- [`get_ambient_cell_map(grid::ManifoldGrid)`]

The following methods have been added to return the invidiual grids
- [`get_topo_grid(grid::ManifoldGrid)`]
- [`get_parametric_grid(grid::ManifoldGrid)`]
- [`get_ambient_grid(grid::ManifoldGrid)`]

The following method has been added to return the panel_ids
- [`get_panel_ids(grid::ManifoldGrid)`]
"""

abstract type ManifoldName end

struct Cube <: ManifoldName end
struct CubedSphere <: ManifoldName end

const cube = Cube()
const cubedsphere = CubedSphere()


struct ManifoldGrid{Dc,Dp,Dp_topo,Dp_parm,A<:Grid{Dc,Dp_topo},B<:Grid{Dc,Dp_parm},C<:Grid{Dc,Dp},E,F,G} <: Grid{Dc,Dp}
  name::ManifoldName
  topo_grid::A
  parametric_grid::B
  ambient_grid::C
  parametric_cell_map::E
  ambient_cell_map::F
  parametric_cell_coords::G
  panel_ids::Vector{Int}
end

function ManifoldGrid(name::ManifoldName,
  topo_grid::Grid{Dc,Dp_topo},
  parametric_grid::Grid{Dc,Dp_parm},
  ambient_grid::Grid{Dc,Dp},
  parametric_cell_map,ambient_cell_map,parametric_cell_coords,panel_ids) where {Dc,Dp,Dp_topo,Dp_parm}

  A = typeof(topo_grid)
  B = typeof(parametric_grid)
  C = typeof(ambient_grid)
  E = typeof(parametric_cell_map)
  F = typeof(ambient_cell_map)
  G = typeof(parametric_cell_coords)

  ManifoldGrid{Dc,Dp,Dp_topo,Dp_parm,A,B,C,E,F,G}(name,topo_grid,parametric_grid,
    ambient_grid,parametric_cell_map,ambient_cell_map,parametric_cell_coords,
    panel_ids)
end

function ManifoldGrid(model::DiscreteModel,name::ManifoldName)
  panel_ids = get_panel_ids(model)
  topo_grid = get_grid(model)
  ManifoldGrid(topo_grid,panel_ids,name)
end

function ManifoldGrid(model::DiscreteModel)
  GenericManifoldGrid(model)
end

function GenericManifoldGrid(model::DiscreteModel)
  println("not implemented yet")
  @notimplemented
end

function cell_maps_from_coords(cell_coords,cell_reffes,cell_type)

  ctype_shapefuns = map(get_shapefuns,cell_reffes)
  cell_shapefuns = expand_cell_data(ctype_shapefuns,cell_type)

  cell_map = lazy_map(linear_combination,cell_coords,cell_shapefuns)

  Fields.MemoArray(cell_map)
end


"""
CubeSurfaceGrid

topo_grid == cubic topology with 3D nodes
parametric_grid == Bump ∘ Rp1 ∘ topo_nodes (only has non-zero nodes on panel 1)
ambient_grid == topo_grid (since topo grid has 3D nodes)
parametric_cell_map == Bump ∘ Rp1 ∘ cmap
ambient_cell_map == cmap

"""
function ManifoldGrid(topo_grid::Grid{Dc,Dp_topo},panel_ids,::Cube) where {Dc,Dp_topo}
  println("cube manifold grid")

  cmaps = get_cell_map(topo_grid)
  cell_node_ids = get_cell_node_ids(topo_grid)
  topo_cell_coords = get_cell_coordinates(topo_grid)

  parametric_cell_coords = get_cube_nodes(topo_cell_coords,panel_ids)
  parametric_nodes = get_panel_1_nodes_from_coords(topo_grid,parametric_cell_coords,panel_ids)

  parametric_grid = Gridap.Geometry.UnstructuredGrid(parametric_nodes,cell_node_ids,
      get_reffes(topo_grid),get_cell_type(topo_grid),OrientationStyle(topo_grid))

  parametric_cell_map = cell_maps_from_coords(parametric_cell_coords, get_reffes(topo_grid),get_cell_type(topo_grid))

  ManifoldGrid(cube,topo_grid,parametric_grid,topo_grid,parametric_cell_map,
              cmaps,parametric_cell_coords,panel_ids)

end


function get_cube_nodes(topo_cell_coords,panel_ids)
  coords_panel1 = lazy_map(PanelMap(), topo_cell_coords, panel_ids)
  coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)
  return coords_panel1_2D
end


"""
CubedSphereGrid

topo_grid == cubic topology with 3D nodes, representing central angles
parametric_grid == Bump ∘ Rp1 ∘ topo_nodes (only has non-zero nodes on panel 1)
ambient_grid == SigmaMap ∘ GnomonicMap ∘ Bump ∘ Rp1 ∘ topo_nodes
parametric_cell_map == Bump ∘ Rp1 ∘ cmap
ambient_cell_map == R1p ∘ SigmaMap ∘ GnomonicMap ∘ Bump ∘ Rp1 ∘ cmap

"""
function ManifoldGrid(topo_grid::Grid{Dc,Dp_topo},panel_ids,::CubedSphere) where {Dc,Dp_topo}
  println("cubed sphere manifold grid")

  cmaps = get_cell_map(topo_grid)
  cell_node_ids = get_cell_node_ids(topo_grid)
  topo_cell_coords = get_cell_coordinates(topo_grid)

  parametric_cell_coords, ambient_cell_coords = get_cubed_sphere_nodes(topo_cell_coords,panel_ids)
  parametric_nodes = get_panel_1_nodes_from_coords(topo_grid,parametric_cell_coords,panel_ids)

  parametric_grid = Gridap.Geometry.UnstructuredGrid(parametric_nodes,cell_node_ids,
      get_reffes(topo_grid),get_cell_type(topo_grid),OrientationStyle(topo_grid))

  parametric_cell_map = cell_maps_from_coords(parametric_cell_coords, get_reffes(topo_grid),get_cell_type(topo_grid))

  ambient_nodes = get_nodes_from_coords(topo_grid,ambient_cell_coords)

  ambient_grid = Gridap.Geometry.UnstructuredGrid(ambient_nodes,cell_node_ids,
      get_reffes(topo_grid),get_cell_type(topo_grid),OrientationStyle(topo_grid))

  ambient_cell_map = cell_maps_from_coords(ambient_cell_coords, get_reffes(topo_grid),get_cell_type(topo_grid))


  ManifoldGrid(cubedsphere,topo_grid,parametric_grid,ambient_grid,parametric_cell_map,
              ambient_cell_map,parametric_cell_coords,panel_ids)

end


function get_cubed_sphere_nodes(topo_cell_coords,panel_ids)

  coords_panel1 = lazy_map(PanelMap(), topo_cell_coords, panel_ids)
  cangles_panel1 = lazy_map(BumpMap(), coords_panel1)

  latlon_panel1 = lazy_map(GnomonicMap(), cangles_panel1)
  sphere_panel1 = lazy_map(Sigma(),latlon_panel1)
  sphere_panelp = lazy_map(InvPanelMap(), sphere_panel1, panel_ids)

  return cangles_panel1, sphere_panelp
end



"""
grid API
"""
Gridap.Geometry.get_cell_map(grid::ManifoldGrid) = grid.parametric_cell_map
Gridap.Geometry.get_cell_coordinates(grid::ManifoldGrid) = grid.parametric_cell_coords
Gridap.Geometry.get_node_coordinates(grid::ManifoldGrid) = get_node_coordinates(grid.parametric_grid)
Gridap.Geometry.get_cell_node_ids(grid::ManifoldGrid) = get_cell_node_ids(grid.parametric_grid)
Gridap.Geometry.get_reffes(grid::ManifoldGrid) = get_reffes(grid.parametric_grid)
Gridap.Geometry.get_cell_type(grid::ManifoldGrid) = get_cell_type(grid.parametric_grid)
Gridap.Geometry.get_cell_ref_coordinates(grid::ManifoldGrid) = get_cell_ref_coordinates(grid.topo_grid)
Gridap.Geometry.OrientationStyle(grid::ManifoldGrid) = Gridap.Geometry.NonOriented()

"""
ambient grid functions
"""
get_ambient_cell_map(grid::ManifoldGrid) = grid.ambient_cell_map
get_ambient_node_coordinates(grid::ManifoldGrid) = get_node_coordinates(grid.ambient_grid)
get_ambient_cell_coordinates(grid::ManifoldGrid) = get_cell_coordinates(grid.ambient_grid)

"""
additional helper functions
"""
get_parametric_grid(grid::ManifoldGrid) = grid.parametric_grid
get_ambient_grid(grid::ManifoldGrid) = grid.ambient_grid
get_topo_grid(grid::ManifoldGrid) = grid.topo_grid

get_panel_ids(grid::ManifoldGrid) = grid.panel_ids
