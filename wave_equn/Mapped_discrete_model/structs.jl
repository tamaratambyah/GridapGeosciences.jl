import Gridap.Geometry.Grid

struct MyMappedGrid{Dc,Dp,A,M,L} <: Grid{Dc,Dp}
  grid::Grid{Dc,Dp}
  geo_map::A
  phys_map::M # New map in the physical space
  node_coords::L
  function MyMappedGrid(grid::Grid{Dc,Dp},phys_map) where {Dc,Dp}

    cell_node_ids = get_cell_node_ids(grid)
    old_nodes = get_node_coordinates(grid)

    dims = num_dims(grid) + 1 # increase dimension by 1
    T = VectorValue{dims,Float64}
    node_coordinates = Vector{T}(undef,length(old_nodes))
    c_coor = get_cell_coordinates(grid)./1

    N = length(c_coor)
    map_c_coor  = fill(Vector{T}(undef,length(c_coor)),N)
    for i = 1:N
      map_c_coor[i] = Base.map( (x)->map(x), c_coor[i])
    end
    Gridap.Geometry._cell_vector_to_dof_vector!(node_coordinates,cell_node_ids,map_c_coor)


    model_map=get_cell_map(grid)
    geo_map= get_cell_map(grid) #lazy_map(∘,phys_map,model_map)
    node_coords = collect(node_coordinates)

    new{Dc,Dp,typeof(geo_map),typeof(phys_map),typeof(node_coords)}(grid,geo_map,phys_map,node_coords)
  end
end


Gridap.Geometry.get_node_coordinates(grid::MyMappedGrid) = grid.node_coords
Gridap.Geometry.get_cell_node_ids(grid::MyMappedGrid) = get_cell_node_ids(grid.grid)
Gridap.Geometry.get_reffes(grid::MyMappedGrid) = get_reffes(grid.grid)
Gridap.Geometry.get_cell_type(grid::MyMappedGrid) = get_cell_type(grid.grid)


struct MyMappedDiscreteModel{Dc,Dp} <: DiscreteModel{Dc,Dp}
  model::DiscreteModel{Dc,Dp}
  mapped_grid
  function MyMappedDiscreteModel(model::DiscreteModel{Dc,Dp},node_coords) where {Dc,Dp}
    new{Dc,Dp}(model,node_coords)
  end
end

Gridap.Geometry.get_grid(model::MyMappedDiscreteModel) = model.mapped_grid
Gridap.Geometry.get_cell_map(model::MyMappedDiscreteModel) = get_cell_map(model.mapped_grid)
Gridap.Geometry.get_grid_topology(model::MyMappedDiscreteModel) = Gridap.Geometry.get_grid_topology(model.model)
Gridap.Geometry.get_face_labeling(model::MyMappedDiscreteModel) = get_face_labeling(model.model)
