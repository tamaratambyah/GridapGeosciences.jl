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
include("cube_topo/cube_surface_1_cell_per_panel.jl")
include("maps/panel_ids_from_refinement_v2.jl")
include("maps/panel_rotations.jl")
include("maps/bump_panel1.jl")
include("maps/maps.jl")


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

cube_model_3D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel()...)

### Coarset model
model = cube_model_3D
# model = Gridap.Adaptivity.refine(model)
panel_ids = get_panel_ids(model)
cell_coords = get_cell_coordinates(model)
cmaps = get_cell_map(model)
cell_node_ids = get_cell_node_ids(model)
ref_cell_coords = get_cell_ref_coordinates(get_grid(model))


### maps required for cube
rotate_panel_p_to_1, rotate_panel_1_to_p = panel_rotations()
rp1 = map(TensorValue,rotate_panel_p_to_1)
r1p = map(TensorValue,rotate_panel_1_to_p)

_A , _B, _b = bump_matrics()
A_bump = TensorValue(_A)
B_bump = TensorValue(_B)
b_bump = VectorValue(_b)


Rp1 = map(p->PanelRotationMap(rp1[p]), 1:6)
R1p = map(p->PanelRotationMap(r1p[p]), 1:6)
Bump = Panel1BumpMap(A_bump,B_bump,b_bump)



struct CellPanelMaps{A,B} <: Map # map from ref FE -> panel p
  Rp1::A
  R1p::A
  Bump::B
end

function Gridap.Arrays.return_cache(f::CellPanelMaps,panel_id::Int64,cmap)
  y = first(f.Rp1)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::CellPanelMaps,panel_id::Int64,cmap)
  y = cache
  y = f.R1p[panel_id] ∘ f.Bump∘ f.Bump∘ f.Rp1[panel_id] ∘ cmap
  # y = Operation(f.R1p[panel_id])(Operation(f.Bump)(Operation(f.Bump)(Operation(f.Rp1[panel_id])(cmap) )  ) )

  # About the same speed. ∘ is easier to read
  return y
end

panel_id = 1
typeof(R1p[panel_id] ∘ Bump∘ Bump∘ Rp1[panel_id] ∘ cmaps[1])


cell_panel_maps = lazy_map(CellPanelMaps(Rp1,R1p,Bump), panel_ids, cmaps)
cache = array_cache(cell_panel_maps)
bm1() = lazy_collect(cache,cell_panel_maps)
@benchmark bm1()

# test_cell_maps(cell_panel_maps,ref_cell_coords,cell_coords)
evaluate(cell_panel_maps[1],ref_cell_coords[1])
# lazy_map(evaluate,cmaps, coord)


# ### this allocates a lot!
# phys_coords = lazy_map(evaluate,cell_panel_maps,ref_cell_coords)
# cache = array_cache(phys_coords)
# bm1() = lazy_collect(cache,phys_coords)
# @benchmark bm1()
