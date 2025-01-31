abstract type ManifoldGrid{Dc,Dp_grid,Dp_topo} <: Grid{Dc,Dp_grid} end

"""
ManifoldGrid -- grid type where Dp_grid can be different from Dp_topo
  Dc      : dimesnion of the reference cells
  Dp_grid : dimesion of the physical grid
  Dp_topo : dimesnion of the topology
"""

"""
CubedSphereGrid --- grid type that contains:
      cube_grid == underlying cube
      sphere_grid == cube_grid mapped to sphere surface
      sphere_cell_map == cells maps from 2D ref FE space → surface of sphere
      cube_to_sphere_map == map from cube → surface of sphere (anlaytical/polynomial)
      transfer_info  == (original analytical map, transfer bool, order)
                        transfer_bool = true, transfer map
                        transfer_bool = false, reinterpolat map

  Dc      : dimesnion of the reference cells  = 2
  Dp_grid : dimesion of the physical grid     = 2 or 3
  Dp_topo : dimesnion of the topology         = 3

  Tp      : Float64                           = eltype(VectorValue{3,Float64})
  O       : typeof( Gridap.Geometry.NonOriented() )
  Tn      : typeof( nothing )
  A       : typeof( sphere_cell_map )
  B       : typeof( cube_to_sphere_map )
"""
struct CubedSphereGrid{Dc,Dp_grid,Dp_topo,Tp,O,Tn,A,B} <: ManifoldGrid{Dc,Dp_grid,Dp_topo} where Tn
  cube_grid::UnstructuredGrid{Dc,Dp_topo,Tp,O,Tn}
  sphere_grid::UnstructuredGrid{Dc,Dp_grid,Tp,O,Tn}
  sphere_cell_map::A
  cube_to_sphere_map::B
  transfer_info::Tuple{Function,Bool,Integer}
end

function CubedSphereGrid(cube_grid::UnstructuredGrid{Dc,Dp_topo,Tp,O,Tn},
  analytical_cube_to_sphere_map::Function) where {Dc,Dp_topo,Tp,O,Tn}

  println("analytical constructor")

  # map the nodes from the cube to the sphere
  cube_nodes = get_node_coordinates(cube_grid)
  sphere_nodes = map( (XYZ)-> analytical_cube_to_sphere_map(XYZ), cube_nodes  ) #lazy_map( map_cube_to_sphere, cube_nodes  )

  sphere_grid = Gridap.Geometry.UnstructuredGrid(sphere_nodes,
                            get_cell_node_ids(cube_grid),
                            get_reffes(cube_grid),
                            get_cell_type(cube_grid),
                            Gridap.Geometry.NonOriented())

  # create map from 2D reference FE to the sphere
  geo_cell_map = Fill( Gridap.Fields.GenericField( analytical_cube_to_sphere_map ), num_cells(cube_grid))
  cube_cell_map = get_cell_map(cube_grid)
  sphere_cell_map = lazy_map(∘,geo_cell_map,cube_cell_map)

  transfer_info = (analytical_cube_to_sphere_map,false,1)

  Dp_grid = num_point_dims(sphere_grid)
  A = typeof(sphere_cell_map)
  B = typeof(analytical_cube_to_sphere_map)
  CubedSphereGrid{Dc,Dp_grid,Dp_topo,Tp,O,Tn,A,B}(cube_grid,sphere_grid,sphere_cell_map,analytical_cube_to_sphere_map,transfer_info)

end

"""
grid API
"""
Gridap.Geometry.OrientationStyle(grid::CubedSphereGrid) = Gridap.Geometry.NonOriented()

function Gridap.Geometry.get_reffes(grid::CubedSphereGrid)
  Gridap.Geometry.get_reffes(grid.sphere_grid)
end

function Gridap.Geometry.get_cell_type(grid::CubedSphereGrid)
  Gridap.Geometry.get_cell_type(grid.sphere_grid)
end

function Gridap.Geometry.get_node_coordinates(grid::CubedSphereGrid)
  Gridap.Geometry.get_node_coordinates(grid.sphere_grid)
end

function Gridap.Geometry.get_cell_node_ids(grid::CubedSphereGrid)
  Gridap.Geometry.get_cell_node_ids(grid.sphere_grid)
end

function Gridap.Geometry.get_cell_map(grid::CubedSphereGrid)
  grid.sphere_cell_map
end

function get_cube_to_sphere_map(grid::CubedSphereGrid)
  grid.cube_to_sphere_map
end

function get_transfer_info(grid::CubedSphereGrid)
  grid.transfer_info
end



abstract type ManifoldDiscreteModel{Dc,Dp_grid,Dp_topo} <:  DiscreteModel{Dc,Dp_grid}  end



struct UnstructuredManifoldDiscreteModel{Dc,Dp_topo,Dp_grid,Tp,O} <: ManifoldDiscreteModel{Dc,Dp_grid,Dp_topo}
  grid::Grid{Dc,Dp_grid}#::UnstructuredGrid{Dc,Dp_grid,Tp,O,Tn}
  grid_topology::UnstructuredGridTopology{Dc,Dp_topo,Tp,O}
  face_labeling::FaceLabeling
end

function ManifoldDiscreteModel(grid::Grid,grid_topology::GridTopology,labels::FaceLabeling)
  UnstructuredManifoldDiscreteModel(grid,grid_topology,labels)
end

# Implementation of the interface

Gridap.Geometry.get_grid(model::UnstructuredManifoldDiscreteModel) = model.grid

Gridap.Geometry.get_grid_topology(model::UnstructuredManifoldDiscreteModel) = model.grid_topology

Gridap.Geometry.get_face_labeling(model::UnstructuredManifoldDiscreteModel) = model.face_labeling



struct AdaptedManifoldDiscreteModel{Dc,Dp_grid,Dp_topo,A<:ManifoldDiscreteModel,
  B<:ManifoldDiscreteModel,C<:AdaptivityGlue} <: ManifoldDiscreteModel{Dc,Dp_grid,Dp_topo}
  model  ::A
  parent ::B
  glue   ::C
end

function AdaptedManifoldDiscreteModel(model,parent::ManifoldDiscreteModel{Dc,Dp_grid,Dp_topo},glue) where {Dc,Dp_grid,Dp_topo}
  @Gridap.Helpers.check !isa(model,AdaptedManifoldDiscreteModel)
  A = typeof(model)
  B = typeof(parent)
  C = typeof(glue)
  return AdaptedManifoldDiscreteModel{Dc,Dp_grid,Dp_topo,A,B,C}(model,parent,glue)
end


# DiscreteModel API
Gridap.Geometry.get_grid(model::AdaptedManifoldDiscreteModel)          = get_grid(model.model)
Gridap.Geometry.get_grid_topology(model::AdaptedManifoldDiscreteModel) = get_grid_topology(model.model)
Gridap.Geometry.get_face_labeling(model::AdaptedManifoldDiscreteModel) = get_face_labeling(model.model)

# Other getters
Gridap.Adaptivity.get_model(model::AdaptedManifoldDiscreteModel)  = model.model
Gridap.Adaptivity.get_parent(model::AdaptedManifoldDiscreteModel{Dc,Dp_grid,Dp_topo,A,<:AdaptedDiscreteModel}) where {Dc,Dp_grid,Dp_topo,A} = get_model(model.parent)
Gridap.Adaptivity.get_parent(model::AdaptedManifoldDiscreteModel{Dc,Dp_grid,Dp_topo,A,B}) where {Dc,Dp_grid,Dp_topo,A,B} = model.parent
Gridap.Adaptivity.get_adaptivity_glue(model::AdaptedManifoldDiscreteModel) = model.glue

# Relationships
"""
Returns true if m1 is a "child" model of m2, i.e., if m1 is the result of adapting m2
"""
function Gridap.Adaptivity.is_child(m1::AdaptedManifoldDiscreteModel,m2::DiscreteModel)
  return get_parent(m1) === m2 # m1 = refine(m2)
end

function Gridap.Adaptivity.is_child(m1::AdaptedManifoldDiscreteModel,m2::AdaptedManifoldDiscreteModel)
  return get_parent(m1) === get_model(m2) # m1 = refine(m2)
end

Gridap.Adaptivity.is_child(m1::DiscreteModel,m2::AdaptedManifoldDiscreteModel) = false
