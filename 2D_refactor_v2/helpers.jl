
function lazy_collect(cache,arr)
  for i in eachindex(arr)
    getindex!(cache, arr, i)
  end
end

function print_lazy_collect(arr,cell_node_ids)
  cache = array_cache(arr)
  for i in eachindex(arr)
    println(getindex!(cache, arr, i))
    println(cell_node_ids[i])
  end
end

function test_coords(arr,cell_coords)
  cache = array_cache(arr)
  for i in eachindex(arr)
    @test getindex!(cache, arr, i) == cell_coords[i]
  end
  println("Yay!")
end



function get_nodes_from_coords(topo::UnstructuredGridTopology{Dc,Dp},
  coords_array::AbstractArray{<:Vector{<:VectorValue{D,T}}}) where {Dc,Dp,D,T}

  cache = array_cache(coords_array)
  cell_node_ids = get_faces(topo,Dc,0)
  nodes = similar(coords_panelp, VectorValue{D,T}, num_vertices(topo))

  for i in eachindex(coords_array)
    ids = cell_node_ids[i]
    nodes[ids] .= getindex!(cache, coords_array, i)
  end

  return nodes
end

function make_grid(topo::UnstructuredGridTopology{Dc,Dp},nodes::AbstractArray) where {Dc,Dp}

  cell_reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  cell_node_ids = get_faces(topo,Dc,0)
  cell_type = get_cell_type(topo)

  Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,OrientationStyle(topo))
end
