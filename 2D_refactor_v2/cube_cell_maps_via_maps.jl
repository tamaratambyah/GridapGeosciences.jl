using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

include("cube_surface_1_cell_per_panel.jl")
include("panel_ids_from_refinement_v2.jl")
include("panel_rotations.jl")
# include("maps.jl")

function lazy_collect(cache,arr)
  for i in eachindex(arr)
    getindex!(cache, arr, i)
  end
end

function _lazy_collect(cache,arr)
  out = similar(arr)
  for i in eachindex(arr)
    out[i] = getindex!(cache, arr, i)
  end
  return out
end

function get_nodes_from_coords(topo::UnstructuredGridTopology{Dc,Dp},coords) where {Dc,Dp}

  cell_node_ids = get_faces(topo,Dc,0)
  nodes = Vector{VectorValue{Dc+1,Float64}}(undef,num_vertices(topo))

  for i in 1:length(cell_node_ids)
    ids = cell_node_ids[i]
    nodes[ids] .= coords[i]
  end

  return nodes
end

function make_grid(topo::UnstructuredGridTopology{Dc,Dp},nodes::AbstractArray) where {Dc,Dp}

  cell_reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  cell_node_ids = get_faces(topo,Dc,0)
  cell_type = get_cell_type(topo)

  Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,OrientationStyle(topo))
end

### mapp_nodes
const rotate_panel_p_to_1, rotate_panel_1_to_p = panel_rotations()

struct PanelMap <: Map # rotate panel p -> panel 1
end

function Gridap.Arrays.return_cache(f::PanelMap,cellx::AbstractArray{<:VectorValue{3}},panel_id::Int64)
  y =similar(cellx) # data from X into y
  A =  [TensorValue( rotate_panel_p_to_1[panel_id] ) for i in 1:length(y)] # rotate panel 1 -> panel p
  return y, A
end

function Gridap.Arrays.evaluate!(cache,f::PanelMap,cellx::AbstractArray{<:VectorValue{3}},panel_id::Int64)
  y, A = cache
  y .= A .⋅ cellx
end


function Gridap.Arrays.return_cache(f::PanelMap,X::VectorValue{3},panel_id::Int64)
  y = VectorValue(get_array(X)) # data from X into y
  A =  TensorValue( rotate_panel_p_to_1[panel_id] ) # rotate panel 1 -> panel p
  return y, A
end

function Gridap.Arrays.evaluate!(cache,f::PanelMap,X::VectorValue{3},panel_id::Int64)
  y, A = cache
  y = A ⋅ X
end



struct InvPanelMap <: Map # rotate panel 1 -> panel p
end

function Gridap.Arrays.return_cache(f::InvPanelMap,cellx::AbstractArray{<:VectorValue{3}},panel_id::Int64)
  y =similar(cellx) # data from X into y
  A =  [TensorValue( rotate_panel_1_to_p[panel_id] ) for i in 1:length(y)] # rotate panel 1 -> panel p
  return y, A
end

function Gridap.Arrays.evaluate!(cache,f::InvPanelMap,cellx::AbstractArray{<:VectorValue{3}},panel_id::Int64)
  y, A = cache
  y .= A .⋅ cellx
end


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

cube_model_3D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel()...)

### Coarset model
model = cube_model_3D
panel_ids = get_panel_ids(model)
cell_coords = get_cell_coordinates(model)
cmaps = get_cell_map(model)

######## Map cell coordinates
coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
cache = array_cache(coords_panel1)
bm1() = lazy_collect(cache,coords_panel1)
@benchmark bm1()

coords_panelp = lazy_map(InvPanelMap(), cell_coords, panel_ids)
cache = array_cache(coords_panelp)
bm2() = lazy_collect(cache,coords_panelp)
@benchmark bm2()

# test apply map and inverse
coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
coords_panelp = lazy_map(InvPanelMap(), coords_panel1, panel_ids)
cache = array_cache(coords_panelp)
bm3() = lazy_collect(cache,coords_panelp)
@benchmark bm3()

# make a new grid
mapp_cell_coords = _lazy_collect(cache,coords_panelp)

mapped_nodes = get_nodes_from_coords(get_grid_topology(model),cell_coords)
mapped_grid = make_grid(get_grid_topology(model),mapped_nodes)
writevtk(mapped_grid,dir*"/cube_model",append=false)

###############################################################################
######## Map cmaps
ref_cell_coords = get_cell_ref_coordinates(get_grid(model))
coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
# Master = lazy_map(Broadcasting(∘),array, cmaps)


_cc = Operation(PanelMap())(cell_coords,panel_ids)

cache = array_cache(Master)
bm2() = lazy_collect(cache,Master)
@benchmark bm2()







### Make another map that takes the cmap as input
struct OtherPanelMap <: Map # rotate panel p -> panel 1
end
function Gridap.Arrays.return_cache(f::OtherPanelMap,cell_ref::AbstractArray{<:VectorValue{2}},panel_id::Int64,cmap)
  cellx = evaluate(cmap,cell_ref)
  y = similar(cellx)
  A =  [TensorValue( rotate_panel_p_to_1[panel_id] ) for i in 1:length(y)] # rotate panel 1 -> panel p
  return y, cellx, A
end

function Gridap.Arrays.evaluate!(cache,f::OtherPanelMap,cell_ref::AbstractArray{<:VectorValue{2}},panel_id::Int64,cmap)
  y, cellx, A = cache
  y .= A .⋅ cellx
end


#
Master = lazy_map(OtherPanelMap(), ref_cell_coords, panel_ids, cmaps)
cache = array_cache(Master)
bm2() = lazy_collect(cache,Master)
@benchmark bm2()
