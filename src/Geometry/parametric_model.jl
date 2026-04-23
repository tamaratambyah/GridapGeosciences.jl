function coarse_parametric_model(radius; a=π/4, npanels=6)
  cube_model = coarse_cube_model(a, npanels)
  panel_model = parametric_model(cube_model,radius)
  return panel_model
end

# TO-DO: decide better name (???)
function setup_panel_cmaps(cube_cmaps, panel_ids)
  cube2panel_map = lazy_map(p->MatMultField( A_cube2panel[p] ), panel_ids)
  panel_cmaps = lazy_map(∘,cube2panel_map,cube_cmaps)
end

function parametric_model(cube_model, radius)
  Dc = num_cell_dims(cube_model)
  Dp = num_point_dims(cube_model)

  @check Dp == Dc+1

  panel_ids = get_panel_ids(cube_model)
  cube_grid = get_grid(cube_model)
  cube_topo = get_grid_topology(cube_model)
  cube_nodes = get_node_coordinates(cube_grid)

  ## map the cube to the parametric domain
  cube_cmaps  = get_cell_map(cube_grid)
  panel_cmaps = setup_panel_cmaps(cube_cmaps, panel_ids)

  ## create the panel model
  ## to correctly trigger Dc=2,Dp=2, need to have 2D nodes.
  ## the panel_nodes are just junk 2D nodes, never used by Gridap
  ## the panel_grid has the bespoke panel_cmap
  ## the panel_topo is the same as the cube, but with the 2D nodes so that Dp=2
  ## the panel_labels are from the panel_topo
  panel_nodes = map(x->Point(x[2],x[3]),cube_nodes) # these are just junk nodes, never used
  panel_grid = Geometry.UnstructuredGrid(panel_nodes,get_cell_node_ids(cube_grid),get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid),
                      nothing,panel_cmaps)
  panel_topo = UnstructuredGridTopology(panel_nodes,get_cell_node_ids(cube_grid),get_cell_type(cube_topo),get_polytopes(cube_topo),OrientationStyle(cube_topo))
  panel_labels = FaceLabeling(panel_topo)

  panel_model = ParametricDiscreteModel(panel_grid,panel_topo,panel_labels,panel_ids,radius)

  return panel_model
end


struct ParametricDiscreteModel{Dc,Dp,Tp,B,Tf} <: DiscreteModel{Dc,Dp}
  grid::UnstructuredGrid{Dc,Dp,Tp,B}
  grid_topology::UnstructuredGridTopology{Dc,Dp,Tp,B}
  face_labeling::FaceLabeling
  panel_ids::AbstractArray{Int}
  forward_map_generator::Tf
end

function ParametricDiscreteModel(grid::UnstructuredGrid{2,2},
                                 grid_topology::UnstructuredGridTopology{2,2},
                                 face_labeling::FaceLabeling,
                                 panel_ids::AbstractArray{Int},
                                 radius::Real)
  forward_map_generator = ForwardMap2DGenerator(radius)
  return ParametricDiscreteModel(grid, grid_topology, face_labeling, panel_ids, forward_map_generator)
end

function ParametricDiscreteModel(grid::UnstructuredGrid{3,3},
                                 grid_topology::UnstructuredGridTopology{3,3},
                                 face_labeling::FaceLabeling,
                                 panel_ids::AbstractArray{Int},
                                 radius::Real,
                                 thickness::Real)
  forward_map_generator = ForwardMap3DGenerator(radius, thickness)
  return ParametricDiscreteModel(grid, grid_topology, face_labeling, panel_ids, forward_map_generator)
end

Geometry.get_grid(model::ParametricDiscreteModel) = model.grid
Geometry.get_grid_topology(model::ParametricDiscreteModel) = model.grid_topology
Geometry.get_face_labeling(model::ParametricDiscreteModel) = model.face_labeling
get_panel_ids(model::ParametricDiscreteModel) = model.panel_ids
Geometry.get_cell_map(model::ParametricDiscreteModel) = model.grid.cell_map
get_forward_map_generator(model::ParametricDiscreteModel) = model.forward_map_generator
get_forward_map_generator(model::AdaptedDiscreteModel{Dc,Dp,<:ParametricDiscreteModel}) where {Dc,Dp} = model.model.forward_map_generator

function get_radius(model::Union{ParametricDiscreteModel,AdaptedDiscreteModel{Dc,Dp,<:ParametricDiscreteModel}}) where {Dc,Dp}
  generator = get_forward_map_generator(model)
  return generator.radius
end

function get_thickness(model::Union{ParametricDiscreteModel{3,3},AdaptedDiscreteModel{3,3,<:ParametricDiscreteModel}})
  generator = get_forward_map_generator(model)
  return generator.thickness
end
