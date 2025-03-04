
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
    @test evaluate(arr[i],ref_coords[i]) == cell_coords[i]
  end
  println("Yay!")
end


function get_nodes_from_coords(topo::UnstructuredGridTopology{Dc,Dp},
  coords_array::AbstractArray{<:Vector{<:VectorValue{D,T}}}) where {Dc,Dp,D,T}

  cache = array_cache(coords_array)
  cell_node_ids = get_faces(topo,Dc,0)
  nodes = similar(coords_array, VectorValue{D,T}, num_vertices(topo))

  for i in eachindex(coords_array)
    ids = cell_node_ids[i]
    nodes[ids] .= getindex!(cache, coords_array, i)
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



function plot_latlons(latlon,simName::String)
  !isdir(plotsdir()) && mkdir(plotsdir())

  markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross]
  _colors = palette(:tab10)
  p1 = plot(title = "Cells")
  p2 = plot(title = "Cell points")


  cache = array_cache(latlon)
  for i in eachindex(latlon)
    panel = panel_ids[i]
    out = getindex!(cache, latlon, i)

    lon = map(x->x[1],out)
    lat = map(x->x[2],out)

    plot!(p1,lon,lat,c=_colors[panel])
    scatter!(p2,lon,lat,marker=markers[panel],c=_colors[panel])
  end
  plot!(p1,legend=false,xlabel="longitude",ylabel="latitude")
  plot!(p2,legend=false,xlabel="longitude",ylabel="latitude")

  savefig(p1,plotsdir()*"/$(simName)_cells")
  savefig(p2,plotsdir()*"/$(simName)_cells_points")
end
