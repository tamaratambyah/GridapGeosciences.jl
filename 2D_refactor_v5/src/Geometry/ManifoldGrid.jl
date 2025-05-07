"""
    ManifoldGrid

Stores the information required to solve on the surface
name                    : name of the manifold (i.e. cube or cubedsphere)
cube_grid_3D            : grid of the associated cube in 3D
parametric_grid         : grid of the parametric domain (for FEM) i.e. panel 1
                          - note, this grid represents panel 1, and only panel 1
                          nodes are non-zero.
ambient_grid            : grid of the ambient space (manifold)
parametric_cell_map     : map to the parametric space (used for integration)
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

The following helper methods have been added to return useful grid information
- [`get_panel_ids(grid::ManifoldGrid)`]
- [`get_manifold_name(grid::ManifoldGrid)`]
"""


abstract type ManifoldName end
struct Cube <: ManifoldName end
struct CubedSphere <: ManifoldName end

const cube = Cube()
const cubedsphere = CubedSphere()


struct ManifoldGrid{Dc,Dp,Dp_amb,A<:Grid{Dc,Dp_amb},B<:Grid{Dc,Dp},C<:Grid{Dc,Dp_amb},E,G} <: Grid{Dc,Dp}
  name::ManifoldName
  cube_grid_3D::A
  parametric_grid::B
  ambient_grid::C
  parametric_cell_map::E
  parametric_cell_coords::G
  panel_ids::Vector{Int}
end


function ManifoldGrid(name::ManifoldName,
  cube_grid_3D::Grid{Dc,Dp_amb},
  parametric_grid::Grid{Dc,Dp},
  ambient_grid::Grid{Dc,Dp_amb},
  parametric_cell_map,parametric_cell_coords,panel_ids) where {Dc,Dp,Dp_amb}

  @check Dp_amb > Dp

  A = typeof(cube_grid_3D)
  B = typeof(parametric_grid)
  C = typeof(ambient_grid)
  E = typeof(parametric_cell_map)
  G = typeof(parametric_cell_coords)

  ManifoldGrid{Dc,Dp,Dp_amb,A,B,C,E,G}(name,cube_grid_3D,parametric_grid,
    ambient_grid,parametric_cell_map,parametric_cell_coords,
    panel_ids)
end

function ManifoldGrid(name::ManifoldName,cube_grid_3D::Grid{Dc,Dp_amb},panel_ids::Vector{Int}) where {Dc,Dp_amb}

  cube_grid_2D, = cube_surface_2D(cube_grid_3D,panel_ids)

  parametric_grid, parametric_cell_map, parametric_cell_coords = construct_parametric_grid(cube_grid_2D,cube_grid_3D,panel_ids)
  ambient_grid = construct_ambient_grid(name,cube_grid_3D,panel_ids)

  ManifoldGrid(name,cube_grid_3D,parametric_grid,ambient_grid,
    parametric_cell_map,parametric_cell_coords,panel_ids)

end


function ManifoldGrid(model::DiscreteModel,name::ManifoldName)
  panel_ids = get_panel_ids(model)

  @assert !all(panel_ids .> 6) """\n
  Invalid panel ids! \n Panel ids must be 1,…,6
  """
  cube_grid_3D = get_grid(model)
  ManifoldGrid(name,cube_grid_3D,panel_ids)
end

function ManifoldGrid(model::DiscreteModel)
  println("not implemented yet")
  @notimplemented
end

"""
parametric_grid

is always a grid relating to [-a,a] on each panel
"""
function construct_parametric_grid(cube_grid_2D::Grid{Dc,Dp},
      cube_grid_3D::Grid{Dc,Dp_amb},panel_ids::Vector{Int}) where {Dc,Dp,Dp_amb}

  cmaps = get_cell_map(cube_grid_3D)

  cell_coords = get_cell_coordinates(cube_grid_3D)
  cell_node_ids = get_cell_node_ids(cube_grid_3D)

  panel1_nodes = get_node_coordinates(cube_grid_2D)


  g =  BumpField(A_bump,B_bump,b_bump)
  k = map(p-> g ∘ PanelRotationField(rp1_3D[p]), panel_ids)
  parametric_cell_map = lazy_map(∘,k,cmaps)

  parametric_cell_coords = get_cube_nodes_2D(cell_coords,panel_ids)

  parametric_grid = UnstructuredGrid(panel1_nodes,cell_node_ids,
      get_reffes(cube_grid_2D),get_cell_type(cube_grid_2D),OrientationStyle(cube_grid_2D),
      nothing,parametric_cell_map)

  return parametric_grid, parametric_cell_map, parametric_cell_coords
end

"""
CubeSurfaceGrid
"""

function construct_ambient_grid(::Cube,cube_grid_3D::Grid{Dc,Dp_amb},panel_ids::Vector{Int}) where {Dc,Dp_amb}
  println("cube manifold grid")

  return cube_grid_3D
end


"""
CubedSphereGrid
"""
function construct_ambient_grid(::CubedSphere,cube_grid_3D::Grid{Dc,Dp_amb},panel_ids::Vector{Int}) where {Dc,Dp_amb}
  println("cubed sphere manifold grid")

  cmaps = get_cell_map(cube_grid_3D)

  cube_cell_coords_3D = get_cell_coordinates(cube_grid_3D)
  ambient_cell_coords = get_cubed_sphere_nodes(cube_cell_coords_3D,panel_ids)

  ambient_nodes = get_nodes_from_coords(cube_grid_3D,ambient_cell_coords)


  g =  BumpField(A_bump,B_bump,b_bump)
  k = map(p-> PanelRotationField(r1p_3D[p]) ∘ SigmaField(r) ∘ GnomonicField() ∘ g ∘ PanelRotationField(rp1_3D[p]), panel_ids)
  ambient_cell_map = lazy_map(∘,k,cmaps)


  ambient_grid = Gridap.Geometry.UnstructuredGrid(ambient_nodes,get_cell_node_ids(cube_grid_3D),
    get_reffes(cube_grid_3D),get_cell_type(cube_grid_3D),OrientationStyle(cube_grid_3D),
    nothing,ambient_cell_map)

  return ambient_grid
end



function get_cubed_sphere_nodes(cube_cell_coords_3D,panel_ids)

  coords_panel1 = lazy_map(Rp1PanelMap3D(), cube_cell_coords_3D, panel_ids)
  cangles_panel1 = lazy_map(BumpMap(), coords_panel1)

  latlon_panel1 = lazy_map(GnomonicMap(), cangles_panel1)
  sphere_panel1 = lazy_map(Sigma(),latlon_panel1)
  sphere_panelp = lazy_map(R1pPanelMap3D(), sphere_panel1, panel_ids)

  return sphere_panelp
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
Gridap.Geometry.get_cell_ref_coordinates(grid::ManifoldGrid) = get_cell_ref_coordinates(grid.cube_grid_3D)
Gridap.Geometry.OrientationStyle(grid::ManifoldGrid) = Gridap.Geometry.NonOriented()


"""
ambient grid functions
"""
get_ambient_node_coordinates(grid::ManifoldGrid) = get_node_coordinates(grid.ambient_grid)
get_ambient_cell_coordinates(grid::ManifoldGrid) = get_cell_coordinates(grid.ambient_grid)
get_ambient_cell_map(grid::ManifoldGrid) = get_cell_map(grid.ambient_grid)


"""
additional helper functions
"""
get_parametric_grid(grid::ManifoldGrid) = grid.parametric_grid
get_ambient_grid(grid::ManifoldGrid) = grid.ambient_grid
get_panel_ids(grid::ManifoldGrid) = grid.panel_ids
get_manifold_name(grid::ManifoldGrid) = grid.name
get_3D_cube_grid(grid::ManifoldGrid) = grid.cube_grid_3D
