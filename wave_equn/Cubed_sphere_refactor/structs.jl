# Dc = num_cell_dims(topo) = dimension of model = 2
# Dp = num_point_dims(grid)  = dimension of points = 3
# Tp = Float64 = VectorValue{3,Float64}
# O = typeof( Gridap.Geometry.NonOriented() )
# A = typeof( sphere_cell_map )
# B = typeof( cube_to_sphere_map )

"""
CubedSphereGrid
"""

struct CubedSphereGrid{Dc,Dp,Tp,O,A,B} <: Grid{Dc,Dp}
  cube_grid::UnstructuredGrid{Dc,Dp,Tp,O}
  sphere_grid::UnstructuredGrid{Dc,Dp,Tp,O}
  sphere_cell_map::A
  cube_to_sphere_map::B
end

"""
CubedSphereGrid -- with analaytical cube_to_sphere_map
"""
function CubedSphereGrid(cube_grid::UnstructuredGrid{Dc,Dp,Tp,O},
    analytical_cube_to_sphere_map::Function) where {Dc,Dp,Tp,O}

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

  A = typeof(sphere_cell_map)
  B = typeof(analytical_cube_to_sphere_map)
  CubedSphereGrid{Dc,Dp,Tp,O,A,B}(cube_grid,sphere_grid,sphere_cell_map,analytical_cube_to_sphere_map)

end

"""
CubedSphereGrid -- with polynomial cube_to_sphere_map (FE_map)
"""
function CubedSphereGrid(cube_grid::UnstructuredGrid{Dc,Dp,Tp,O},
  analytical_cube_to_sphere_map::Function, order::Integer ) where {Dc,Dp,Tp,O}

  map_basis = get_fe_basis(get_fe_space(FE_map))
  map_trian = get_triangulation(map_basis)
  map_reffes = get_reffes(map_trian)
  order = get_order(map_reffes[1])

  cube_model = UnstructuredDiscreteModel(cube_grid)

  # make a scalar FE
  V_scal = FESpace(cube_model,
                  ReferenceFE(lagrangian,Float64,order);conformity=:H1)

  cell_node_ids = get_cell_dof_ids(V_scal)
  sphere_nodes = Vector{VectorValue{3,Float64}}(undef,num_free_dofs(V_scal))

  c_dofs = get_fe_dof_basis(V_scal)
  ref_nodes = lazy_map(get_nodes,get_data(c_dofs))
  vhx = lazy_map(evaluate,get_data(FE_map),ref_nodes)


  Gridap.Geometry._cell_vector_to_dof_vector!(sphere_nodes,cell_node_ids,vhx)


  sphere_grid = Gridap.Geometry.UnstructuredGrid(sphere_nodes,
                            cell_node_ids,
                            get_reffes(cube_grid),
                            get_cell_type(cube_grid),
                            Gridap.Geometry.NonOriented())

  sphere_cell_map = get_cell_map(sphere_grid)

  A = typeof(sphere_cell_map)
  B = typeof(FE_map)
  CubedSphereGrid{Dc,Dp,Tp,O,A,B}(cube_grid,sphere_grid,sphere_cell_map,FE_map)

end



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



"""
CSDiscreteModel
"""
struct CubedSphereDiscreteModel{Dc,Dp,Tp,O,A,B} <: DiscreteModel{Dc,Dp}
  CSgrid::CubedSphereGrid{Dc,Dp,Tp,O,A,B}
  grid_topology::UnstructuredGridTopology{Dc,Dp,Tp,O}
  face_labeling::FaceLabeling
end

function Gridap.Geometry.get_grid(model::CubedSphereDiscreteModel)
  model.CSgrid
end


function Gridap.Geometry.get_grid_topology(model::CubedSphereDiscreteModel)
  model.grid_topology
end


function Gridap.Geometry.get_face_labeling(model::CubedSphereDiscreteModel)
  model.face_labeling
end
