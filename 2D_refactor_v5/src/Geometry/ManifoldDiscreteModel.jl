"""
    ManifoldDiscreteModel

Stores the ManifoldGrid in a model structure

manifold_grid     : bespoke manifold grid  (i.e. cube or cubed sphere)
grid_topology     : topology of the underlying cube
face_labeling     : face labeling
parametric_model  : model associated to parametric grid
ambient_model     : model associated to ambient grid

The model interface is overloaded to return
- [`get_grid(model::ManifoldDiscreteModel)`]
- [`get_grid_topology(model::ManifoldDiscreteModel)`]
- [`get_face_labeling(model::ManifoldDiscreteModel)`]

 The following method has been added to return useful infomation
- [`get_panel_ids(model::ManifoldDiscreteModel)`]
- [`get_manifold_name(model::ManifoldDiscreteModel)`]
- [`get_parametric_model(model::ManifoldDiscreteModel)`]
- [`get_ambient_model(model::ManifoldDiscreteModel)`]
"""

struct ManifoldDiscreteModel{Dc,Dp,Dp_amb,A<:Grid{Dc,Dp},B<:GridTopology{Dc,Dp},
          E<:DiscreteModel{Dc,Dp},F<:DiscreteModel{Dc,Dp_amb}} <: DiscreteModel{Dc,Dp}
  manifold_grid:: A
  grid_topology:: B
  face_labeling:: FaceLabeling
  parametric_model:: E
  ambient_model:: F
end

function ManifoldDiscreteModel(manifold_grid::ManifoldGrid{Dc,Dp},topo::GridTopology{Dc,Dp},
  labels::FaceLabeling,parametric_model::DiscreteModel{Dc,Dp},ambient_model::DiscreteModel{Dc,Dp_amb}) where {Dc,Dp,Dp_amb}
  A = typeof(manifold_grid)
  B = typeof(topo)
  E = typeof(parametric_model)
  F = typeof(ambient_model)
  ManifoldDiscreteModel{Dc,Dp,Dp_amb,A,B,E,F}(manifold_grid,topo,labels,parametric_model,ambient_model)
end

function ManifoldDiscreteModel(manifold_grid::ManifoldGrid{Dc,Dp},topo_2D::GridTopology{Dc,Dp},
  labels::FaceLabeling) where {Dc,Dp}

  parametric_model = UnstructuredDiscreteModel(get_parametric_grid(manifold_grid),topo_2D,labels)
  ambient_grid = get_ambient_grid(manifold_grid)
  ambient_topo = UnstructuredGridTopology(ambient_grid)
  ambient_model = UnstructuredDiscreteModel(ambient_grid,ambient_topo,labels)

  ManifoldDiscreteModel(manifold_grid,topo_2D,labels,parametric_model,ambient_model)
end

function ManifoldDiscreteModel(manifold_grid::ManifoldGrid{Dc,Dp}) where {Dc,Dp}

  panel_ids = get_panel_ids(manifold_grid)
  cube_grid_3D = get_3D_cube_grid(manifold_grid)
  Dp_amb = num_point_dims(cube_grid_3D)
  @check Dp_amb > Dp

  cube_grid_2D,topo_2D,labels = cube_surface_2D(cube_grid_3D,panel_ids)

  parametric_model = UnstructuredDiscreteModel(get_parametric_grid(manifold_grid),topo_2D,labels)
  ambient_grid = get_ambient_grid(manifold_grid)
  ambient_topo = UnstructuredGridTopology(ambient_grid)
  ambient_model = UnstructuredDiscreteModel(ambient_grid,ambient_topo,labels)

  ManifoldDiscreteModel(manifold_grid,topo_2D,labels,parametric_model,ambient_model)
end


function ManifoldDiscreteModel(model::DiscreteModel,name::ManifoldName)
  manifold_grid = ManifoldGrid(model,name)
  ManifoldDiscreteModel(manifold_grid)
end

Gridap.Geometry.get_grid(model::ManifoldDiscreteModel) = model.manifold_grid
Gridap.Geometry.get_grid_topology(model::ManifoldDiscreteModel) = model.grid_topology
Gridap.Geometry.get_face_labeling(model::ManifoldDiscreteModel) = model.face_labeling
Gridap.Geometry.get_cell_coordinates(model::ManifoldDiscreteModel) = get_cell_coordinates(get_grid(model))

get_panel_ids(model::ManifoldDiscreteModel) = get_panel_ids(get_grid(model))
get_manifold_name(model::ManifoldDiscreteModel) = get_manifold_name(get_grid(model))
get_parametric_model(model::ManifoldDiscreteModel) = model.parametric_model
get_ambient_model(model::ManifoldDiscreteModel) = model.ambient_model


get_manifold_name(amodel::AdaptedDiscreteModel) = get_manifold_name(amodel.model)
get_parametric_model(amodel::AdaptedDiscreteModel) = get_parametric_model(amodel.model)
get_ambient_model(amodel::AdaptedDiscreteModel) = get_ambient_model(amodel.model)



function get_latlon_model(model::ManifoldDiscreteModel)
  @check get_manifold_name(model) == CubedSphere()
  ambient_model = get_ambient_model(model)
  ambient_grid = get_grid(ambient_model)

  ambient_nodes = get_node_coordinates(ambient_grid)

  latlon_nodes = collect1d(lazy_map(Sigma(),ambient_nodes))

  latlon_grid = UnstructuredGrid(latlon_nodes,get_cell_node_ids(ambient_grid),
  get_reffes(ambient_grid),get_cell_type(ambient_grid),OrientationStyle(ambient_grid))

  UnstructuredDiscreteModel(latlon_grid,get_grid_topology(model),get_face_labeling(model))

end

get_latlon_model(amodel::AdaptedDiscreteModel) = get_latlon_model(amodel.model)
