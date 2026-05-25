# AtlasGrids.jl
#
# Defines AtlasGrid{Dc,Da}: a Grid{Dc,Dc} whose cells carry local reference coords
# (Dc-dimensional, one coordinate system per coarse chart).  Physical Da-dimensional
# coords are NOT stored here; they are computed only in visualization_data.

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Helpers
using Gridap.Adaptivity, Gridap.Visualization
using FillArrays

# ============================================================
# AtlasGrid
# ============================================================

"""
    AtlasGrid{Dc,Da,G,A,M,O} <: Gridap.Geometry.Grid{Dc,Dc}

A DG-style grid for a `Dc`-dimensional manifold atlas.  Each fine cell stores its
corners in the **local reference frame** of its coarse chart (Dc-dimensional).
Physical Da-dimensional coordinates are never materialised here; they are computed
on demand only during visualization.

# Fields
- `param_grid`        — fine `Grid{Dc,Dc}` from uniform refinement (topology/connectivity).
- `cell_local_coords` — `Table{Point{Dc}}`: per-cell local corner coords, DG-style.
- `cell_to_chart`     — per-fine-cell index of the originating coarse chart (1-based).
- `physical_maps`     — one callable per coarse chart: `Point{Dc} → Point{Da}`.
- `orientation_style` — kept explicitly because atlas orientation can differ from the
                        underlying grid topology (e.g. a Möbius strip).
"""
struct AtlasGrid{Dc, Da,
                 G <: Gridap.Geometry.Grid{Dc,Dc},
                 A <: AbstractVector,
                 M,
                 O <: Gridap.Geometry.OrientationStyle} <: Gridap.Geometry.Grid{Dc,Dc}
  param_grid        :: G
  cell_local_coords :: A
  cell_to_chart     :: AbstractVector{<:Integer}
  physical_maps     :: M
  orientation_style :: O

  function AtlasGrid(
    param_grid        :: Gridap.Geometry.Grid{Dc,Dc},
    cell_local_coords :: AbstractVector,
    cell_to_chart     :: AbstractVector{<:Integer},
    physical_maps,
    orientation_style :: Gridap.Geometry.OrientationStyle,
  ) where Dc
    # Infer Da from the first physical map evaluation (cache pattern: one call, constructor only)
    sample_pt  = cell_local_coords[1][1]
    fwd0       = physical_maps[cell_to_chart[1]]
    cache0     = Gridap.Arrays.return_cache(fwd0, sample_pt)
    sample_out = Gridap.Arrays.evaluate!(cache0, fwd0, sample_pt)
    Da = length(sample_out)
    n = Gridap.Geometry.num_cells(param_grid)
    @check length(cell_local_coords) == n
    @check length(cell_to_chart) == n
    G = typeof(param_grid)
    A = typeof(cell_local_coords)
    M = typeof(physical_maps)
    O = typeof(orientation_style)
    new{Dc,Da,G,A,M,O}(param_grid, cell_local_coords,
                        cell_to_chart, physical_maps, orientation_style)
  end
end

# ----------------------------------------------------------
# Outer constructor: coarse_model + coarse_local_coords + physical_maps
# ----------------------------------------------------------

"""
    _build_atlas_grid(coarse_model, coarse_local_coords, physical_maps,
                      num_refinements, orientation_style) -> (AtlasGrid, fine_model)

Private helper shared by the `AtlasGrid` and `AtlasDiscreteModel` outer constructors.
Performs the uniform refinement exactly once and returns both the `AtlasGrid` and the
`fine_model` so the caller can extract topology and face labeling without re-refining.
"""
function _build_atlas_grid(
    coarse_model        :: Gridap.Geometry.DiscreteModel{Dc,Dc},
    coarse_local_coords,
    physical_maps,
    num_refinements :: Int,
    orientation_style,
) where Dc

  # ── 1. Refine coarse_model N times → fine_model, cell_to_chart ───────────
  current_model = coarse_model
  cell_to_chart = collect(1:num_cells(coarse_model))
  first_glue    = nothing

  for l in 1:num_refinements
    adapted       = Gridap.Adaptivity.refine(current_model)
    level_map     = adapted.glue.n2o_faces_map[Dc+1]
    cell_to_chart = cell_to_chart[level_map]
    if l == 1
      first_glue = adapted.glue
    end
    current_model = adapted.model
  end
  fine_model = current_model

  # ── 2. Get N-level reference positions from Gridap's internal ref-grid ───
  ref_model = Gridap.Adaptivity.get_ref_grid(first_glue.refinement_rules[1])
  for _ in 2:num_refinements
    ref_model = Gridap.Adaptivity.refine(ref_model).model
  end
  ref_coords  = Gridap.Geometry.get_cell_coordinates(Gridap.Geometry.get_grid(ref_model))
  n_per_chart = Gridap.Geometry.num_cells(Gridap.Geometry.get_grid(ref_model))

  # ── 3. Build per-chart local maps Ψ_k : ref_element → coarse_local_coords[k] ──
  reffe     = Gridap.Geometry.get_reffes(Gridap.Geometry.get_grid(coarse_model))[1]
  shapefuns = Gridap.ReferenceFEs.get_shapefuns(reffe)
  Ψ_maps    = map(corners -> linear_combination(corners, shapefuns), coarse_local_coords)

  # ── 4. Tile ref_coords and apply Ψ_k — one cache per cell, reused per corner ──
  ncharts      = num_cells(coarse_model)
  child_ids    = repeat(1:n_per_chart, ncharts)
  ref_per_fine = lazy_map(Reindex(ref_coords), child_ids)
  chart_Ψ      = lazy_map(Reindex(Ψ_maps), cell_to_chart)
  local_lazy   = lazy_map(chart_Ψ, ref_per_fine) do Ψ, corners
    cache = Gridap.Arrays.return_cache(Ψ, corners[1])
    map(c -> Gridap.Arrays.evaluate!(cache, Ψ, c), corners)
  end
  cell_local_coords_table = Gridap.Arrays.Table(collect.(local_lazy))

  # ── 5. Resolve orientation_style ─────────────────────────────────────────
  os = if isnothing(orientation_style)
    Gridap.Geometry.OrientationStyle(Gridap.Geometry.get_grid(fine_model))
  else
    orientation_style
  end

  atlas_grid = AtlasGrid(
    Gridap.Geometry.get_grid(fine_model),
    cell_local_coords_table,
    cell_to_chart,
    physical_maps,
    os,
  )

  atlas_grid, fine_model
end

"""
    AtlasGrid(coarse_model, coarse_local_coords, physical_maps, num_refinements=1;
              orientation_style=nothing)

Build an `AtlasGrid` by:
1. Refining `coarse_model` `num_refinements` times (for topology + cell_to_chart).
2. Obtaining N-level reference positions from Gridap's internal `get_ref_grid` machinery
   (general: works for any polytope).
3. Building per-chart local maps `Ψ_k` from `coarse_local_coords[k]` using isoparametric
   interpolation (`linear_combination` with shape functions) — no auxiliary model needed.
4. Assembling `cell_local_coords` via `lazy_map` with explicit cache reuse (one cache per
   cell, reused over all corners of that cell).

`coarse_local_coords[k]` gives the local corner coordinates of coarse cell k in whatever
local frame the user defines for that chart.  Different charts may have different local
coordinate systems.

If `orientation_style` is `nothing`, it is copied from the fine param_grid.
"""
function AtlasGrid(
    coarse_model        :: Gridap.Geometry.DiscreteModel{Dc,Dc},
    coarse_local_coords,
    physical_maps,
    num_refinements :: Int = 1;
    orientation_style = nothing,
) where Dc
  atlas_grid, _ = _build_atlas_grid(
    coarse_model, coarse_local_coords, physical_maps, num_refinements, orientation_style)
  atlas_grid
end

# ----------------------------------------------------------
# Gridap.Geometry.Grid{Dc,Dc} interface
# ----------------------------------------------------------

Gridap.Geometry.OrientationStyle(
  ::Type{<:AtlasGrid{Dc,Da,G,A,M,O}}) where {Dc,Da,G,A,M,O} = O()

Gridap.Geometry.num_cells(g::AtlasGrid) = length(g.cell_to_chart)

Gridap.Geometry.get_reffes(g::AtlasGrid) = Gridap.Geometry.get_reffes(g.param_grid)

Gridap.Geometry.get_cell_type(g::AtlasGrid) = Gridap.Geometry.get_cell_type(g.param_grid)

Gridap.Geometry.get_cell_coordinates(g::AtlasGrid) = g.cell_local_coords

Gridap.Geometry.get_node_coordinates(g::AtlasGrid) = g.cell_local_coords.data

function Gridap.Geometry.get_cell_node_ids(g::AtlasGrid)
  cc = g.cell_local_coords
  Gridap.Arrays.Table(Int32.(1:length(cc.data)), cc.ptrs)
end

# ----------------------------------------------------------
# Custom API
# ----------------------------------------------------------

get_param_grid(g::AtlasGrid)                      = g.param_grid
get_physical_maps(g::AtlasGrid)                   = g.physical_maps
get_cell_to_chart(g::AtlasGrid)                   = g.cell_to_chart
get_ambient_dim(::AtlasGrid{Dc,Da}) where {Dc,Da} = Da

# ============================================================
# Coarse mesh library and convenience constructors
# ============================================================

include("CoarseMeshes.jl")

"""
    AtlasGrid(info::CoarseMeshInfo, num_refinements; orientation_style=nothing)

Build an `AtlasGrid` from a `CoarseMeshInfo`, using the physical maps stored in `info`.
"""
function AtlasGrid(
    info            :: CoarseMeshInfo,
    num_refinements :: Int;
    orientation_style = nothing,
)
  AtlasGrid(info.model, info.local_coords, info.physical_maps, num_refinements;
            orientation_style)
end

"""
    AtlasGrid(info::CoarseMeshInfo, custom_maps, num_refinements; orientation_style=nothing)

Build an `AtlasGrid` from a `CoarseMeshInfo`, overriding the default physical maps.
"""
function AtlasGrid(
    info            :: CoarseMeshInfo,
    custom_maps,
    num_refinements :: Int;
    orientation_style = nothing,
)
  AtlasGrid(info.model, info.local_coords, custom_maps, num_refinements;
            orientation_style)
end

"""
    AtlasGrid(mesh::AbstractCoarseMesh, num_refinements; orientation_style=nothing)

Build an `AtlasGrid` directly from a mesh descriptor (e.g. `CubedSphereMesh(1.0)`).
Calls `get_coarse_mesh(mesh)` internally and uses the default physical maps.
"""
function AtlasGrid(
    mesh            :: AbstractCoarseMesh,
    num_refinements :: Int;
    orientation_style = nothing,
)
  AtlasGrid(get_coarse_mesh(mesh), num_refinements; orientation_style)
end
