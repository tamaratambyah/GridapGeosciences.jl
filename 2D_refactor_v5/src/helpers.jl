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

function test_cell_maps(model::DiscreteModel)
  test_cell_maps(get_cell_map(model),
              get_cell_ref_coordinates(get_grid(model)),
              get_cell_coordinates(get_grid(model)))

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
  coords_array::AbstractArray,
  panel_ids::Vector{Int}) where {Dc,Dp}

  T = eltype(first(testitem(coords_array)))
  D = length(first(testitem(coords_array)))

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



function plot_coords(coords,panel_ids,cell_node_ids;plotTitle::String="test",xlab::String="x",ylab::String="y")
  !isdir(plotsdir()) && mkdir(plotsdir())

  markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross]
  _colors = palette(:tab10)

  ids = [1,2,4,3] # for plotting Gridap ordered nodes
  plot()

  cache = array_cache(coords)
  for i in eachindex(coords)
    panel = panel_ids[i]
    out = getindex!(cache, coords, i)
    q = text.((cell_node_ids[i])[ids], halign=:left, valign=:bottom,10)
    plot!(extract_VectorValue(out)...,seriestype=:path,linestyle=:solid,lw=2,
    series_annotations = q,
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


function true_area(name::ManifoldName)
  if name == CubedSphere()
    4*π*RADIUS^2
  elseif name == Cube()
    6*(2*1.0)^2
  end
end

function get_refined_models(model,levs_refinement::Int)
  models = Any[model]
  for i in 1:levs_refinement
    ref_model = Adaptivity.refine(models[i])
    push!(models,ref_model)
  end
  models
end


function get_surface_area(manifold_model,degree::Int)
  Ω = Triangulation(manifold_model)
  m = Metric(manifold_model)
  dΩg = Measure(m,Ω,degree)

  computed_area = sum( integrate(1.0,dΩg))

  return computed_area
end


l2(e,dΩ) = sum(∫(e⊙e)dΩ)



"""
Round a cellwise array of vector values
"""
struct RoundVectorValues <: Map
end

function Gridap.Arrays.return_cache(f::RoundVectorValues,
  cellx::AbstractArray{<:VectorValue})
  x = first(cellx)
  T = typeof(x)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::RoundVectorValues,
  cellx::AbstractArray{<:VectorValue})
  y = cache
  map!(x -> VectorValue(round.(get_array(x))) , y, cellx)
  return y
end




function my_atan(x,y)

  if x == 0.0 && y == 0.0
    println("centre")
    return angle = 0.0
  end


  # x axis
  if y == 0.0
    println("x axis")
    if x > 0.0
      return angle = 0.0
    elseif x < 0.0
      return angle = π
    end
  end

  # y axis
  if x == 0.0
    println("y axis")
    if y > 0.0
      return angle = π/2
    elseif y < 0.0
      return angle = 3*π/2
    end
  end

  ref_angle = atan(abs(y)/abs(x))
  if x > 0.0 && y > 0.0 # quad 1
    println("quad 1")
    angle = ref_angle
  elseif x < 0.0 && y > 0.0 # quad 2
    println("quad 2")
    angle = π -  ref_angle
  elseif x < 0.0 && y < 0.0 # quad 3
    println("quad 3")
    angle = π +  ref_angle
  elseif x > 0.0 && y < 0.0 # quad 4
    println("quad 4")
    angle = 2*π -  ref_angle
  end

  return angle
end




struct ExtractVectorValues{A} <: Map
  id::A
end

function Gridap.Arrays.return_cache(f::ExtractVectorValues,
  cellx::AbstractArray{<:VectorValue})
  x = first(cellx)
  T = typeof(x[f.id])
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::ExtractVectorValues,
  cellx::AbstractArray{<:VectorValue})
  y = cache
  map!(x -> x[f.id] , y, cellx)
  return y
end
