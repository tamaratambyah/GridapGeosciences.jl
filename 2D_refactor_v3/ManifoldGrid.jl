"""
    ManifoldGrid

Stores the information required to solve on the surface

topo_grid         : grid of the topology (i.e. cube)
phys_grid         : grid of the physical domain (for FEM) i.e. panel 1
                    - note, this grid represents panel 1, and only panel 1 nodes
                    are non-zero.
ambient_grid      : grid of the ambient space (manifold)
phys_cell_map     : map to the physical space
ambient_cell_map  : map to ambient space (manifold - not the integration space)
phys_cell_coords  : the cell coords in physical space
                    - this must be stored separatately as we cannot have dG type
                    nodes
panel_ids         : panel ids corresponding to the panels of the cube

The grid interface is overloaded to return the ambient_grid
- [`get_node_coordinates(grid::ManifoldGrid)]
- [`get_cell_node_ids(grid::ManifoldGrid)`]
- [`get_reffes(grid::ManifoldGrid)`]
- [`get_cell_type(grid::ManifoldGrid)`]
- [`get_cell_ref_coordinates(grid::ManifoldGrid)`]
- [`get_cell_coordinates(grid::ManifoldGrid)`]
- ['OrientationStyle(grid::ManifoldGrid)``]

The following methods have been added to return the phys_grid
- [`get_phys_grid(grid::ManifoldGrid)`]
- [`get_phys_node_coordinates(grid::ManifoldGrid)`]
- [`get_phys_cell_coordinates(grid::ManifoldGrid)`]
- [`get_phys_cell_map(grid::ManifoldGrid)`]

The following methods have been added to return the topo_grid
- [`get_topo_grid(grid::ManifoldGrid)`]

The following methods have been added to return the panel_ids
- [`get_panel_ids(grid::ManifoldGrid)`]
"""

struct ManifoldGrid{Dc,Dp,Dp_topo,Dp_phys,A<:Grid{Dc,Dp_topo},B<:Grid{Dc,Dp_phys},C<:Grid{Dc,Dp},E,F,G,H} <: Grid{Dc,Dp}
  topo_grid::A
  phys_grid::B
  ambient_grid::C
  phys_cell_map::E
  ambient_cell_map::F
  phys_cell_coords::G
  panel_ids::H
end

function ManifoldGrid(topo_grid::Grid{Dc,Dp_topo},phys_grid::Grid{Dc,Dp_phys},
  ambient_grid::Grid{Dc,Dp},phys_cell_map,ambient_cell_map,phys_cell_coords,
  panel_ids) where {Dc,Dp,Dp_topo,Dp_phys}

  A = typeof(topo_grid)
  B = typeof(phys_grid)
  C = typeof(ambient_grid)
  E = typeof(phys_cell_map)
  F = typeof(ambient_cell_map)
  G = typeof(phys_cell_coords)
  H = typeof(panel_ids)

  ManifoldGrid{Dc,Dp,Dp_topo,Dp_phys,A,B,C,E,F,G,H}(topo_grid,phys_grid,ambient_grid,
    phys_cell_map,ambient_cell_map,phys_cell_coords,panel_ids)
end

function ManifoldGrid(model::DiscreteModel)
  topo_grid = get_grid(model)
  panel_ids = get_panel_ids(model)
  GenericManifoldGrid()
end

function GenericManifoldGrid()
  println("not implemented yet")
end


"""
CubeSurfaceGrid

topo_grid == cubic topology with 3D nodes
phys_grid == Bump ∘ Rp1 ∘ topo_nodes (only has non-zero nodes on panel 1)
ambient_grid == topo_grid (since topo grid has 3D nodes)
phys_cell_map == Bump ∘ Rp1 ∘ cmap
ambient_cell_map == cmap

"""
function CubeGrid(topo_grid::Grid{Dc,Dp_topo},panel_ids) where {Dc,Dp_topo}
  cmaps = get_cell_map(topo_grid)
  cell_node_ids = get_cell_node_ids(topo_grid)
  topo_cell_coords = get_cell_coordinates(topo_grid)

  phys_cell_coords = get_cube_nodes(topo_cell_coords,panel_ids)
  phys_nodes = get_panel_1_nodes_from_coords(topo_grid,phys_cell_coords,panel_ids)

  phys_grid = Gridap.Geometry.UnstructuredGrid(phys_nodes,cell_node_ids,
      get_reffes(topo_grid),get_cell_type(topo_grid),OrientationStyle(topo_grid))

  phys_cell_map = lazy_map(CubePhysCellMap(), panel_ids, cmaps)

  ManifoldGrid(topo_grid,phys_grid,topo_grid,phys_cell_map,cmaps,phys_cell_coords,panel_ids)

end

function CubeGrid(model::DiscreteModel)
  panel_ids = get_panel_ids(model)
  topo_grid = get_grid(model)
  CubeGrid(topo_grid,panel_ids)
end

function get_cube_nodes(topo_cell_coords,panel_ids)
  coords_panel1 = lazy_map(PanelMap(), topo_cell_coords, panel_ids)
  coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)
  # coords_panel1_3D = lazy_map(BumpMap(), coords_panel1_2D)
  # coords_panelp = lazy_map(InvPanelMap(), coords_panel1_3D, panel_ids)
  return coords_panel1_2D
end


"""
CubedSphereGrid

topo_grid == cubic topology with 3D nodes, representing central angles
phys_grid == Bump ∘ Rp1 ∘ topo_nodes (only has non-zero nodes on panel 1)
ambient_grid == SigmaMap ∘ GnomonicMap ∘ Bump ∘ Rp1 ∘ topo_nodes
phys_cell_map == Bump ∘ Rp1 ∘ cmap
ambient_cell_map == R1p ∘ SigmaMap ∘ GnomonicMap ∘ Bump ∘ Rp1 ∘ cmap

"""
function CubedSphereGrid(topo_grid::Grid{Dc,Dp_topo},panel_ids) where {Dc,Dp_topo}
  cmaps = get_cell_map(topo_grid)
  cell_node_ids = get_cell_node_ids(topo_grid)
  topo_cell_coords = get_cell_coordinates(topo_grid)

  phys_cell_coords, ambient_cell_coords = get_cubed_sphere_nodes(topo_cell_coords,panel_ids)
  phys_nodes = get_panel_1_nodes_from_coords(topo_grid,phys_cell_coords,panel_ids)

  phys_grid = Gridap.Geometry.UnstructuredGrid(phys_nodes,cell_node_ids,
      get_reffes(topo_grid),get_cell_type(topo_grid),OrientationStyle(topo_grid))

  phys_cell_map = lazy_map(CubePhysCellMap(), panel_ids, cmaps)



  ambient_nodes = get_nodes_from_coords(topo_grid,ambient_cell_coords)

  ambient_grid = Gridap.Geometry.UnstructuredGrid(ambient_nodes,cell_node_ids,
      get_reffes(topo_grid),get_cell_type(topo_grid),OrientationStyle(topo_grid))

  ambient_cell_map = lazy_map(SphereAmbientCellMap(), panel_ids, cmaps)


  ManifoldGrid(topo_grid,phys_grid,ambient_grid,phys_cell_map,ambient_cell_map,phys_cell_coords,panel_ids)

end

function CubedSphereGrid(model::DiscreteModel)
  panel_ids = get_panel_ids(model)
  topo_grid = get_grid(model)
  CubedSphereGrid(topo_grid,panel_ids)
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
function Gridap.Geometry.get_node_coordinates(grid::ManifoldGrid)
  Gridap.Geometry.get_node_coordinates(grid.ambient_grid)
end

function get_phys_node_coordinates(grid::ManifoldGrid)
  Gridap.Geometry.get_node_coordinates(grid.phys_grid)
end

function Gridap.Geometry.get_cell_node_ids(grid::ManifoldGrid)
  Gridap.Geometry.get_cell_node_ids(grid.ambient_grid)
end

function Gridap.Geometry.get_reffes(grid::ManifoldGrid)
  Gridap.Geometry.get_reffes(grid.ambient_grid)
end

function Gridap.Geometry.get_cell_type(grid::ManifoldGrid)
  Gridap.Geometry.get_cell_type(grid.ambient_grid)
end

function Gridap.Geometry.get_cell_map(grid::ManifoldGrid)
  grid.ambient_cell_map
end

function get_phys_cell_map(grid::ManifoldGrid)
  grid.phys_cell_map
end

function get_phys_grid(grid::ManifoldGrid)
  grid.phys_grid
end

function get_topo_grid(grid::ManifoldGrid)
  grid.topo_grid
end

function Gridap.Geometry.get_cell_ref_coordinates(grid::ManifoldGrid)
  get_cell_ref_coordinates(grid.topo_grid)
end

function Gridap.Geometry.get_cell_coordinates(grid::ManifoldGrid)
  get_cell_coordinates(grid.ambient_grid)
end

function get_phys_cell_coordinates(grid::ManifoldGrid)
  grid.phys_cell_coords
end

Gridap.Geometry.OrientationStyle(grid::ManifoldGrid) = Gridap.Geometry.NonOriented()

function get_panel_ids(grid::ManifoldGrid)
  grid.panel_ids
end
