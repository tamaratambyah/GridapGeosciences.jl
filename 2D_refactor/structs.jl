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

Gridap.Geometry.num_point_dims(model::UnstructuredManifoldDiscreteModel) = num_point_dims(get_grid(model))

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
struct CubedSphereGrid{Dc,Dp_grid,Dp_topo,Tp,O,Tn,A,B} <: Grid{Dc,Dp_grid} where Tn
  cube_grid::UnstructuredGrid{Dc,Dp_topo,Tp,O,Tn}
  sphere_grid::UnstructuredGrid{Dc,Dp_grid,Tp,O,Tn}
  sphere_cell_map::A
  cube_to_sphere_map::B
  transfer_info::Tuple{Function,Bool,Integer}
end

"""
CubedSphereGrid -- with analaytical cube_to_sphere_map
"""
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
CubedSphereGrid -- convert analytical map to polynomial map
"""
function CubedSphereGrid(cube_grid::UnstructuredGrid{Dc,Dp_topo,Tp,O,Tn},
  analytical_cube_to_sphere_map::Function,order::Integer; transfer::Bool=false ) where {Dc,Dp_topo,Tp,O,Tn}

  cube_model = UnstructuredDiscreteModel(cube_grid)
  T_vec = eltype(get_node_coordinates(cube_model))
  V_vec = FESpace(cube_model,
                  ReferenceFE(lagrangian,T_vec,order),
                  conformity=:H1)
  FE_map = interpolate(analytical_cube_to_sphere_map,V_vec)

  sphere_grid = get_sphere_grid_polynomial_mapping(cube_model,cube_grid,FE_map,order)
  sphere_cell_map = get_cell_map(sphere_grid)
  transfer_info = (analytical_cube_to_sphere_map,transfer,order)

  Dp_grid = num_point_dims(sphere_grid)
  A = typeof(sphere_cell_map)
  B = typeof(FE_map)
  CubedSphereGrid{Dc,Dp_grid,Dp_topo,Tp,O,Tn,A,B}(cube_grid,sphere_grid,sphere_cell_map,FE_map,transfer_info)

end


"""
CubedSphereGrid -- polynomial map on refined model
"""
function CubedSphereGrid(cube_modelh::AdaptedManifoldDiscreteModel,maph::FEFunction,order::Integer,transfer_info::Tuple{Function,Bool,Integer})
  println("refined CS with polynomial map")

  cube_gridh =  get_grid(cube_modelh)
  Dc = num_cell_dims(cube_grid)
  Dp_topo = num_point_dims(cube_grid)
  Tp = eltype(eltype(get_node_coordinates(cube_grid)))
  O = typeof(OrientationStyle(cube_grid))
  Tn = typeof(nothing)

  sphere_grid = get_sphere_grid_polynomial_mapping(cube_modelh,cube_gridh,maph,order)
  sphere_cell_map = get_cell_map(sphere_grid)

  Dp_grid = num_point_dims(sphere_grid)
  A = typeof(sphere_cell_map)
  B = typeof(maph)
  CubedSphereGrid{Dc,Dp_grid,Dp_topo,Tp,O,Tn,A,B}(cube_grid,sphere_grid,sphere_cell_map,maph,transfer_info)

end


"""
get_sphere_grid_polynomial_mapping -- applyies polynomial mapping to cube nodes
"""
function get_sphere_grid_polynomial_mapping(cube_model,cube_grid,FE_map::FEFunction,order::Integer)

  T_vec = eltype(get_node_coordinates(cube_model)) # VectorValue{3,Float64}

  # make a scalar FE
  T_scal = eltype(T_vec) # Float64
  V_scal = FESpace(cube_model,
                  ReferenceFE(lagrangian,T_scal,order);conformity=:H1)

  cell_node_ids = get_cell_dof_ids(V_scal)
  sphere_nodes = Vector{T_vec}(undef,num_free_dofs(V_scal))

  c_dofs = get_fe_dof_basis(V_scal)
  ref_nodes = lazy_map(get_nodes,get_data(c_dofs))
  vhx = lazy_map(evaluate,get_data(FE_map),ref_nodes)


  Gridap.Geometry._cell_vector_to_dof_vector!(sphere_nodes,cell_node_ids,vhx)


  sphere_grid = Gridap.Geometry.UnstructuredGrid(sphere_nodes,
                            cell_node_ids,
                            get_reffes(cube_grid),
                            get_cell_type(cube_grid),
                            Gridap.Geometry.NonOriented())
  return sphere_grid
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
