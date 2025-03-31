function lazy_collect(arr)
  cache = array_cache(arr)
  lazy_collect(cache,arr)
end

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

function print_lazy_arr(arr)
  cache = array_cache(arr)
  for i in eachindex(arr)
    println(getindex!(cache, arr, i))
  end
end

function test_coords(arr,cell_coords)
  cache = array_cache(arr)
  for i in eachindex(arr)
    @test getindex!(cache, arr, i) == cell_coords[i]
  end
  println("Yay!")
end


function test_cell_maps(arr,ref_coords,cell_coords)
  for i in eachindex(arr)
    # @test (evaluate(arr[i],ref_coords[i]) .- cell_coords[i]) .< 1e-12
    @test norm(evaluate(arr[i],ref_coords[i]) - cell_coords[i], Inf) < 1e-12
  end
  println("Yay!")
end



function get_nodes_from_coords(topo::UnstructuredGridTopology{Dc,Dp},
  coords_array::AbstractArray{<:Vector{<:VectorValue{D,T}}}) where {Dc,Dp,D,T}

  cell_node_ids = get_faces(topo,Dc,0)
  nodes = similar(coords_array, VectorValue{D,T}, num_vertices(topo))

  get_nodes_from_coords!(nodes,cell_node_ids,coords_array)

  return nodes
end

function get_nodes_from_coords(grid::Grid{Dc,Dp},
  coords_array::AbstractArray{<:Vector{<:VectorValue{D,T}}}) where {Dc,Dp,D,T}

  cell_node_ids = get_cell_node_ids(grid)
  nodes = similar(coords_array, VectorValue{D,T}, num_nodes(grid))

  get_nodes_from_coords!(nodes,cell_node_ids,coords_array)

  return nodes
end

function get_nodes_from_coords!(nodes,cell_node_ids,
  coords_array::AbstractArray{<:Vector{<:VectorValue}})

  cache = array_cache(coords_array)

  for i in eachindex(coords_array)
    ids = cell_node_ids[i]
    nodes[ids] .= getindex!(cache, coords_array, i)
  end

end

function get_panel_1_nodes_from_coords(grid::Grid{Dc,Dp},
  coords_array::AbstractArray{<:Vector{<:VectorValue{D,T}}},
  panel_ids::Vector{Int}) where {Dc,Dp,D,T}

  # nodes = similar(coords_array, VectorValue{D,T}, num_nodes(grid))
  nodes = zeros(VectorValue{D,T}, num_nodes(grid))
  cell_node_ids = get_cell_node_ids(grid)

  cache = array_cache(coords_array)

  for i in eachindex(coords_array)
    ids = cell_node_ids[i]
    if panel_ids[i] == 1
      nodes[ids] .= getindex!(cache, coords_array, i)
    end
  end

  return nodes

end


function make_grid(topo::UnstructuredGridTopology{Dc,Dp},coords_array::AbstractArray) where {Dc,Dp}
  nodes = get_nodes_from_coords(topo,coords_array)

  cell_reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  cell_node_ids = get_faces(topo,Dc,0)
  cell_type = get_cell_type(topo)

  Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,OrientationStyle(topo))
end



function plot_coords(coords,panel_ids;plotTitle::String="test",xlab::String="x",ylab::String="y")
  !isdir(plotsdir()) && mkdir(plotsdir())

  markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross]
  _colors = palette(:tab10)

  # ids = [1,2,4,3,1] # for plotting Gridap ordered nodes
  plot()

  cache = array_cache(coords)
  for i in eachindex(coords)
    panel = panel_ids[i]
    out = getindex!(cache, coords, i)

    plot!(extract_VectorValue(out)...,seriestype=:path,linestyle=:solid,lw=2,
          c=_colors[panel],marker=markers[panel])
  end

  plot!(legend=false,xlabel="$xlab",ylabel="$ylab",title="$plotTitle")
  savefig(plotsdir()*"/$(plotTitle)_mesh")
end


function extract_VectorValue(out::AbstractVector{VectorValue{D,Float64}}) where {D}
  ids = [1,2,4,3,1] # for plotting Gridap ordered nodes
  if D == 2
    x = map(x->x[1],out)
    y = map(x->x[2],out)
    return x[ids], y[ids]
  elseif D == 3
    x = map(x->x[1],out)
    y = map(x->x[2],out)
    z = map(x->x[3],out)
    return x[ids], y[ids], z[ids]
  end
end
