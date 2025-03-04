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
include("helpers.jl")
include("cube_surface_1_cell_per_panel.jl")
include("panel_ids_from_refinement_v2.jl")
include("panel_rotations.jl")
include("bump_panel1.jl")


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

cube_model_3D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel()...)

### Coarset model
model = cube_model_3D
panel_ids = get_panel_ids(model)
cell_coords = get_cell_coordinates(model)
cmaps = get_cell_map(model)
cell_node_ids = get_cell_node_ids(model)



################################################################################
##### Rotation map
struct PanelRotationMap{A} <: Map # rotation panel p using mats[p]
  mats::A
end

function Gridap.Arrays.return_cache(f::PanelRotationMap,cellx,panel_id::Int64)
  A = f.mats[panel_id]
  x = first(cellx)
  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationMap,cellx,panel_id::Int64)
  y = cache
  A = f.mats[panel_id]
  map!(x -> A⋅x, y, cellx)
  return y
end

rotate_panel_p_to_1, rotate_panel_1_to_p = panel_rotations()
rp1 = map(TensorValue,rotate_panel_p_to_1)
r1p = map(TensorValue,rotate_panel_1_to_p)

# rotate panel p -> panel 1
PanelMap() = PanelRotationMap(rp1)

coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
cache = array_cache(coords_panel1)
bm1() = lazy_collect(cache,coords_panel1)
@benchmark bm1()

print_lazy_collect(coords_panel1,cell_node_ids)

# rotate panel 1 -> panel p
InvPanelMap() = PanelRotationMap(r1p)

coords_panelp = lazy_map(InvPanelMap(), coords_panel1, panel_ids)
cache = array_cache(coords_panelp)
bm2() = lazy_collect(cache,coords_panelp)
@benchmark bm2()

print_lazy_collect(coords_panelp,cell_node_ids)


################################################################################
##### Bump map
struct Panel1BumpMap{A,B,b} <: Map
  A_bump::A
  B_bump::B
  b_bump::b
end

function Gridap.Arrays.return_cache(f::Panel1BumpMap,cellx::AbstractArray{<:VectorValue{D}}) where {D}
  A = f.A_bump
  B = f.B_bump
  x = first(cellx)

  if D == 3 # D==3, -> y == 2 components; bump 3D -> 2D
    T = typeof(A⋅x)
  elseif D == 2 # D==2 -> y == 3 components;  bump 2D -> 3D
    T = typeof(B⋅x)
  end

  y = similar(cellx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Panel1BumpMap,cellx::AbstractArray{<:VectorValue{D}}) where {D}
  y = cache
  A = f.A_bump
  B = f.B_bump
  b = f.b_bump

  if D == 3 # bump 3D -> 2D
    # y .= A.⋅cellx
    map!(x -> A⋅x, y, cellx)
  elseif D == 2 # bump 2D -> 3D
    # y .= B.⋅cellx .+ b
    map!(x -> B⋅x+b, y, cellx)
  end

  return y
end

_A , _B, _b = bump_matrics()
A_bump = TensorValue(_A)
B_bump = TensorValue(_B)
b_bump = VectorValue(_b)

BumpMap() = Panel1BumpMap(A_bump,B_bump,b_bump)

# take panel p 3D coords, rotate to panel 1, bump to 2D, bump to 3D, rotate to panel p
# should recover original cell_coords
coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)
coords_panel1_3D = lazy_map(BumpMap(), coords_panel1_2D)
coords_panelp = lazy_map(InvPanelMap(), coords_panel1_3D, panel_ids)

cache = array_cache(coords_panelp)
bm() = lazy_collect(cache,coords_panelp)
@benchmark bm()

test_coords(coords_panelp,cell_coords)

print_lazy_collect(coords_panelp,cell_node_ids)
print_lazy_collect(coords_panel1_2D,cell_node_ids)


# ################################################################################
# make a new grid from mapped points

mapped_nodes = get_nodes_from_coords(get_grid_topology(model),coords_panelp)
mapped_grid = make_grid(get_grid_topology(model),mapped_nodes)
writevtk(mapped_grid,dir*"/cube_model_3D",append=false)


# mapped_nodes = get_nodes_from_coords(get_grid_topology(model),coords_panel1_2D)
# mapped_grid = make_grid(get_grid_topology(model),mapped_nodes)
# writevtk(mapped_grid,dir*"/cube_model_ref_2D",append=false)

###############################################################################
#### test on refined model

model = Gridap.Adaptivity.refine(cube_model_3D)
panel_ids = get_panel_ids(model)
cell_coords = get_cell_coordinates(model)
cmaps = get_cell_map(model)
cell_node_ids = get_cell_node_ids(model)


coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)
coords_panel1_3D = lazy_map(BumpMap(), coords_panel1_2D)
coords_panelp = lazy_map(InvPanelMap(), coords_panel1_3D, panel_ids)

cache = array_cache(coords_panelp)
bm() = lazy_collect(cache,coords_panelp)
@benchmark bm()

test_coords(coords_panelp,cell_coords)


print_lazy_collect(coords_panelp,cell_node_ids)

mapped_nodes = get_nodes_from_coords(get_grid_topology(model),coords_panelp)
mapped_grid = make_grid(get_grid_topology(model),mapped_nodes)
writevtk(mapped_grid,dir*"/ref_cube_model_3D",append=false)


# ###############################################################################
# ######## Map cmaps
# ref_cell_coords = get_cell_ref_coordinates(get_grid(model))


# cell_coords = lazy_map(evaluate,cmaps,ref_cell_coords)
# panel1_coords_3D = lazy_map(PanelMap(), cell_coords, panel_ids)



# Master = lazy_map(Broadcasting(∘),array, cmaps)
# cache = array_cache(Master)
# bm2() = lazy_collect(cache,Master)
# @benchmark bm2()

# evaluate(Master[1],ref_cell_coords[1])
# # Broadcasting(Operation(∘))(array,cmaps)
