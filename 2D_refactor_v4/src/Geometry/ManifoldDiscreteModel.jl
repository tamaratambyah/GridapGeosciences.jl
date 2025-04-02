"""
    ManifoldDiscreteModel

Stores the ManifoldGrid in a model structure

manifold_grid : bespoke manifold grid  (i.e. cube or cubed sphere)
grid_topology : topology of the underlying cube
face_labeling : face labeling

The model interface is overloaded to return
- [`get_grid(model::ManifoldDiscreteModel)`]
- [`get_grid_topology(model::ManifoldDiscreteModel)`]
- [`get_face_labeling(model::ManifoldDiscreteModel)`]

 The following method has been added to return useful infomation
- [`get_panel_ids(model::ManifoldDiscreteModel)`]
- [`get_manifold_name(model::ManifoldDiscreteModel)`]
"""


struct ManifoldDiscreteModel{Dc,Dp,Dp_topo,A<:Grid{Dc,Dp},B<:GridTopology{Dc,Dp_topo}} <: DiscreteModel{Dc,Dp}
  manifold_grid:: A
  grid_topology:: B
  face_labeling:: FaceLabeling
end

function ManifoldDiscreteModel(manifold_grid::ManifoldGrid,topo::GridTopology{Dc,Dp_topo},
          labels::FaceLabeling) where {Dc,Dp_topo}
  A = typeof(manifold_grid)
  B = typeof(topo)
  Dp = num_point_dims(manifold_grid)
  ManifoldDiscreteModel{Dc,Dp,Dp_topo,A,B}(manifold_grid,topo,labels)
end

function ManifoldDiscreteModel(manifold_grid::ManifoldGrid)
  topo_grid = get_topo_grid(manifold_grid)
  topo = UnstructuredGridTopology(topo_grid)
  face_labels = FaceLabeling(topo)
  ManifoldDiscreteModel(manifold_grid,topo,face_labels)
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
get_manifold_name(model::AdaptedDiscreteModel) = get_manifold_name(model.model)

"""
helpers to make DiscreteModel for the ambient space
"""
function AmbientDiscreteModel(model::ManifoldDiscreteModel)
  manifold_grid = get_grid(model)
  ambient_grid = get_ambient_grid(manifold_grid)
  topo = get_grid_topology(manifold_model)
  labels = get_face_labeling(manifold_model)
  Geometry.GenericDiscreteModel(ambient_grid,topo,labels)
end

function AmbientDiscreteModel(amodel::AdaptedDiscreteModel)
  model = amodel.model
  AmbientDiscreteModel(model)
end
