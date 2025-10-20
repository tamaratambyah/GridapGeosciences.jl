using Gridap
using GridapMakie
using GLMakie
using GridapGeosciences
using DrWatson

domain = (0, 1, 0, 1)
cell_nums = (10, 10)
model = CartesianDiscreteModel(domain, cell_nums)
Ω = Triangulation(model)

fig = plot(Ω)
wireframe!(Ω, color=:black, linewidth=2)
scatter!(Ω, marker=:star8, markersize=20, color=:blue)


include("Laplace/analytic_funcs.jl")
include("convergence_tools.jl")

s_models = get_refined_models(3)

model = s_models[2]
Ω_panel = Triangulation(model)
panel_ids = get_panel_ids(model)
uh = panelwise_cellfield(f_sin,Ω_panel,panel_ids)

plot(uh)
fig, _ , plt = plot(Ω_panel, uh)


function latlon_geo_map_func(panel_ids::AbstractArray{Int})
  println("latolon serial geo map")
  return lazy_map(p -> Cartesian2SphereicalMap() ∘ MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
end

using Gridap.Geometry
topo = Gridap.Geometry.get_grid_topology(model)
grid = get_grid(model)

nodes =  Gridap.Geometry.get_node_coordinates(grid)

## map the cube to the parametric domain
cmaps  = get_cell_map(grid)
ambient = latlon_geo_map_func(panel_ids)
ambient_cmaps = lazy_map(∘,ambient,cmaps)
pts = get_cell_ref_coordinates(grid)

function get_nodes_from_coords(grid::Grid{Dc,Dp},
  coords_array::AbstractArray{<:Vector{<:VectorValue{D,T}}}) where {Dc,Dp,D,T}

  cell_node_ids = get_cell_node_ids(grid)
  nodes = similar(coords_array, VectorValue{D,T}, Gridap.Geometry.num_nodes(grid))

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

coords = lazy_map(evaluate,ambient_cmaps,pts)
_nodes = get_nodes_from_coords(grid,coords)

_grid = Geometry.UnstructuredGrid(_nodes,get_cell_node_ids(grid),get_reffes(grid),get_cell_type(grid),OrientationStyle(grid),
                    nothing,ambient_cmaps)
_topo = UnstructuredGridTopology(_nodes,get_cell_node_ids(grid),get_cell_type(topo),get_polytopes(topo),OrientationStyle(topo))
_labels = FaceLabeling(_topo)

latlon_model = UnstructuredDiscreteModel(_grid,_topo,_labels)

Ω = Triangulation(latlon_model)
plot(Ω)
wireframe!(Ω, color=:black, linewidth=2)
