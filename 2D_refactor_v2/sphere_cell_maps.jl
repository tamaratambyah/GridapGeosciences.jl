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

include("initialise.jl")
include("sphere_nodes_via_panel1.jl")

### Coarset model
_model = ref_ref_model
panel_ids = get_panel_ids(_model)
cell_coords,  = cubed_sphere_nodes(_model)

cmaps = get_cell_map(_model)
cell_node_ids = get_cell_node_ids(_model)
ref_cell_coords = get_cell_ref_coordinates(get_grid(_model))


### maps required for cube
Rp1 = map(p->PanelRotationMap(rp1[p]), 1:6)
R1p = map(p->PanelRotationMap(r1p[p]), 1:6)
Bump = Panel1BumpMap(A_bump,B_bump,b_bump)
Cangels = CentralAngleMap()


struct SphereCellPanelMaps{A,B,C} <: Map # map from ref FE -> panel p
  Rp1::A
  R1p::A
  Bump::B
  Cangels::C
end

function Gridap.Arrays.return_cache(f::SphereCellPanelMaps,panel_id::Int64,cmap)
  y = first(f.Rp1)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SphereCellPanelMaps,panel_id::Int64,cmap)
  y = cache
  y =  f.Cangels ∘ f.Bump ∘ f.Rp1[panel_id] ∘ cmap

  # About the same speed. ∘ is easier to read
  return y
end



cell_panel_maps = lazy_map(SphereCellPanelMaps(Rp1,R1p,Bump,Cangels), panel_ids, cmaps)
cache = array_cache(cell_panel_maps)
bm1() = lazy_collect(cache,cell_panel_maps)
@benchmark bm1()

### because panel 1 is unifrom in cangels, this test fails

test_cell_maps(cell_panel_maps,ref_cell_coords,cell_coords)

# lazy_map(evaluate,cmaps, coord)


# ### this allocates a lot!
# phys_coords = lazy_map(evaluate,cell_panel_maps,ref_cell_coords)
# cache = array_cache(phys_coords)
# bm1() = lazy_collect(cache,phys_coords)
# @benchmark bm1()
