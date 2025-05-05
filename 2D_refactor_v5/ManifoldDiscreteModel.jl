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

function ManifoldDiscreteModel(model::DiscreteModel,name::ManifoldName)
  manifold_grid = ManifoldGrid(model,name)
  topo = get_grid_topology(model)
  labels = get_face_labeling(model)

  parametric_model = UnstructuredDiscreteModel(get_parametric_grid(manifold_grid),topo,labels)

  ambient_grid = get_ambient_grid(manifold_grid)
  ambient_topo = UnstructuredGridTopology(ambient_grid)
  ambient_model = UnstructuredDiscreteModel(ambient_grid,ambient_topo,labels)

  ManifoldDiscreteModel(manifold_grid,topo,labels,parametric_model,ambient_model)
end

Gridap.Geometry.get_grid(model::ManifoldDiscreteModel) = model.manifold_grid
Gridap.Geometry.get_grid_topology(model::ManifoldDiscreteModel) = model.grid_topology
Gridap.Geometry.get_face_labeling(model::ManifoldDiscreteModel) = model.face_labeling
Gridap.Geometry.get_cell_coordinates(model::ManifoldDiscreteModel) = get_cell_coordinates(get_grid(model))

get_panel_ids(model::ManifoldDiscreteModel) = get_panel_ids(get_grid(model))
get_manifold_name(model::ManifoldDiscreteModel) = get_manifold_name(get_grid(model))
get_manifold_name(model::AdaptedDiscreteModel) = get_manifold_name(model.model)
get_parametric_model(model::ManifoldDiscreteModel) = model.parametric_model
get_ambient_model(model::ManifoldDiscreteModel) = model.ambient_model
