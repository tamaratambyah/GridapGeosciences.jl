# AtlasGrids.jl
#
# Defines AtlasGrid{Dc,Da}: a Grid{Dc,Dc} whose cells carry local reference coords
# (Dc-dimensional, one coordinate system per coarse chart).  Physical Da-dimensional
# coords are NOT stored here; they are computed only in visualization_data.

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Helpers
using Gridap.Adaptivity, Gridap.Visualization
import Gridap.TensorValues: symmetric_part

# ============================================================
# IdentityField
# ============================================================

"""
    IdentityField{D}()

Zero-cost `Field` that maps every `Point{D}` to itself.
Used as the physical map when no embedding is available (intrinsic manifold
without an explicit physical-space representation).

`gradient(IdentityField{D}())` returns `one(TensorValue{D,D,Float64})` at every
point — the Jacobian of the identity is the identity tensor.
"""
struct IdentityField{D} <: Field end

Gridap.Arrays.evaluate!(cache, ::IdentityField{D}, x::Point{D,T}) where {D,T} = x

function Gridap.Arrays.return_cache(::IdentityField{D}, xs::AbstractArray{<:Point{D,T}}) where {D,T}
  similar(xs)
end
function Gridap.Arrays.evaluate!(cache, ::IdentityField{D}, xs::AbstractArray{<:Point{D}}) where D
  copyto!(cache, xs)
  cache
end

function Gridap.Arrays.evaluate!(
    cache, ::FieldGradient{1,IdentityField{D}}, x::Point{D}) where D
  one(TensorValue{D,D,Float64})
end
function Gridap.Arrays.return_cache(
    ::FieldGradient{1,IdentityField{D}}, xs::AbstractArray{<:Point{D}}) where D
  CachedArray(similar(xs, TensorValue{D,D,Float64}))
end
function Gridap.Arrays.evaluate!(
    cache, ::FieldGradient{1,IdentityField{D}}, xs::AbstractArray{<:Point{D}}) where D
  setsize!(cache, size(xs))
  fill!(cache.array, one(TensorValue{D,D,Float64}))
  cache.array
end

# ============================================================
# ManifoldStyle trait
# ============================================================

"""
    ManifoldStyle

Trait controlling how `get_cell_map` assembles the cell map in `AtlasGrid`.

- `ExtrinsicManifold()` — compose RefFE→Chart with Chart→Physical (default).
  Gridap integration uses the full 3D Jacobian; the metric is implicit.
- `IntrinsicManifold()` — return only the RefFE→Chart map.
  The metric must be queried explicitly via `get_cell_metric` to write PDEs.
"""
abstract type ManifoldStyle end
struct ExtrinsicManifold <: ManifoldStyle end
struct IntrinsicManifold <: ManifoldStyle end

# ============================================================
# JtJ: metric from Jacobian transpose  (g = JᵀJ, J = ∇φ)
# ============================================================

"""
    JtJ <: Map

`Map` that computes the pullback metric `g = symmetric_part(J ⋅ Jᵀ)` from the
Jacobian transpose `J = ∇φ` (a `TensorValue{Dc,Da}`).  Used with `Operation` to
build per-cell metric fields lazily:

    metric_field = Operation(JtJ())(gradient(physical_map))
"""
struct JtJ <: Gridap.Arrays.Map end
Gridap.Arrays.evaluate!(cache, ::JtJ, J) = symmetric_part(J ⋅ transpose(J))

# ============================================================
# AtlasGrid
# ============================================================

"""
    AtlasGrid{Dc,Da,G,A,P,C,O,M} <: Gridap.Geometry.Grid{Dc,Dc}

A DG-style grid for a `Dc`-dimensional manifold atlas.  Each fine cell stores its
corners in the **local reference frame** of its coarse chart (Dc-dimensional).
Physical Da-dimensional coordinates are never materialised here; they are computed
on demand only during visualization.

# Type parameters
- `Dc` — cell dimension (= manifold dimension)
- `Da` — ambient (physical) space dimension
- `G`  — concrete type of `param_grid`
- `A`  — concrete type of `cell_local_coords`
- `P`  — concrete type of `cell_physical_maps`
- `C`  — concrete type of `cell_metric`
- `O`  — concrete `OrientationStyle` subtype
- `M`  — concrete `ManifoldStyle` subtype (drives `get_cell_map` dispatch)

# Fields
- `param_grid`          — fine `Grid{Dc,Dc}` from uniform refinement (topology/connectivity).
- `cell_local_coords`   — per-cell local corner coords, DG-style.  May be a lazy
                          `LazyArray` (serial path via `_build_atlas_grid`) or an
                          eager `Table` (distributed path via p4est).
- `cell_physical_maps`  — lazy per-cell map `Point{Dc} → Point{Da}`.  Always concrete:
                          `Fill(IdentityField{Dc}(), n)` when no embedding is provided,
                          giving `Da = Dc` and chart-space visualization.
- `cell_metric`         — lazy per-cell `Field{Dc → SymTensorValue{Dc}}`.  Always concrete:
                          either an analytic expression from `CoarseMeshInfo` or
                          `Operation(JtJ())(gradient(φ))` computed lazily from the
                          physical map.
- `orientation_style`   — kept explicitly because atlas orientation can differ from the
                          underlying grid topology (e.g. a Möbius strip).
- (M)                   — `ManifoldStyle` encoded as type parameter; not stored as a field.
                          Access via `ManifoldStyle(g)`.
"""
struct AtlasGrid{Dc, Da,
                 G <: Gridap.Geometry.Grid{Dc,Dc},
                 A <: AbstractVector,
                 P <: AbstractVector,
                 C <: AbstractVector,
                 O <: Gridap.Geometry.OrientationStyle,
                 M <: ManifoldStyle} <: Gridap.Geometry.Grid{Dc,Dc}
  param_grid         :: G
  cell_local_coords  :: A
  cell_physical_maps :: P
  cell_metric        :: C
  orientation_style  :: O

  function AtlasGrid(
    param_grid         :: Gridap.Geometry.Grid{Dc,Dc},
    cell_local_coords  :: AbstractVector,
    cell_physical_maps :: AbstractVector,
    cell_metric        :: AbstractVector,
    orientation_style  :: Gridap.Geometry.OrientationStyle,
    manifold_style     :: ManifoldStyle,
  ) where Dc
    sample_pt = cell_local_coords[1][1]
    fwd0      = cell_physical_maps[1]
    Da        = length(zero(Gridap.Arrays.return_type(fwd0, sample_pt)))
    n = Gridap.Geometry.num_cells(param_grid)
    @check length(cell_local_coords)  == n
    @check length(cell_physical_maps) == n
    @check length(cell_metric)        == n
    G = typeof(param_grid)
    A = typeof(cell_local_coords)
    P = typeof(cell_physical_maps)
    C = typeof(cell_metric)
    O = typeof(orientation_style)
    M = typeof(manifold_style)
    new{Dc,Da,G,A,P,C,O,M}(param_grid, cell_local_coords, cell_physical_maps, cell_metric, orientation_style)
  end
end

# ----------------------------------------------------------
# ManifoldStyle accessor (mirrors OrientationStyle pattern)
# ----------------------------------------------------------

ManifoldStyle(::Type{<:AtlasGrid{Dc,Da,G,A,P,C,O,M}}) where {Dc,Da,G,A,P,C,O,M} = M()
ManifoldStyle(g::AtlasGrid) = ManifoldStyle(typeof(g))

# ----------------------------------------------------------
# Private helper: build RefFE→Chart maps (common to both styles)
# ----------------------------------------------------------

function _chart_maps(g::AtlasGrid)
  reffes         = Gridap.Geometry.get_reffes(g)
  cell_types     = Gridap.Geometry.get_cell_type(g)
  shapefuns      = map(Gridap.ReferenceFEs.get_shapefuns, reffes)
  cell_shapefuns = lazy_map(Reindex(shapefuns), cell_types)
  lazy_map(linear_combination, g.cell_local_coords, cell_shapefuns)
end

# ----------------------------------------------------------
# _build_atlas_grid: shared private constructor
# ----------------------------------------------------------

"""
    _build_atlas_grid(coarse_model, coarse_local_coords, physical_maps, metric_fields,
                      num_refinements, orientation_style, manifold_style) -> (AtlasGrid, fine_model)

Private helper shared by the `AtlasGrid` and `AtlasDiscreteModel` outer constructors.
Refines `coarse_model` `num_refinements` times and returns both the `AtlasGrid` and the
`fine_model` so the caller can extract topology and face labeling without re-refining.

`num_refinements` must be ≥ 1.  The face labeling of `coarse_model` is propagated to
`fine_model` by Gridap's refinement machinery.

# Arguments
- `coarse_model`       — coarse `DiscreteModel{Dc,Dc}` with face labeling.
- `coarse_local_coords`— `AbstractVector`; entry `k` is a vector of `Point{Dc}` giving
                         the corners of chart `k` in its local reference frame.
- `physical_maps`      — `AbstractVector` of `Field`s; `physical_maps[k]` maps
                         `Point{Dc} → Point{Da}` for chart `k`.
- `metric_fields`      — `AbstractVector` of `Field`s; `metric_fields[k]` maps
                         `Point{Dc} → SymTensorValue{Dc}` (the pullback metric for chart `k`).
- `num_refinements`    — number of uniform refinements (≥ 1).
- `orientation_style`  — grid orientation; `nothing` → copied from the fine param_grid.
- `manifold_style`     — `ExtrinsicManifold()` or `IntrinsicManifold()`.
"""
function _build_atlas_grid(
    coarse_model        :: Gridap.Geometry.DiscreteModel{Dc,Dc},
    coarse_local_coords,
    physical_maps,
    metric_fields,
    num_refinements :: Int,
    orientation_style,
    manifold_style,
) where Dc

  @check num_refinements >= 1 "num_refinements must be ≥ 1, got $num_refinements"

  # ── 1. Refine coarse_model N times → fine_model, cell_to_chart ───────────
  ncharts       = Gridap.Geometry.num_cells(coarse_model)
  cell_to_chart = collect(1:ncharts)

  adapted       = Gridap.Adaptivity.refine(coarse_model)
  first_glue    = adapted.glue
  cell_to_chart = cell_to_chart[first_glue.n2o_faces_map[Dc+1]]
  current_model = adapted.model

  for _ in 2:num_refinements
    adapted       = Gridap.Adaptivity.refine(current_model)
    cell_to_chart = cell_to_chart[adapted.glue.n2o_faces_map[Dc+1]]
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

  # ── 4. Tile ref_coords and apply Ψ_k lazily ──────────────────────────────────
  child_ids     = repeat(1:n_per_chart, ncharts)
  ref_per_fine  = lazy_map(Reindex(ref_coords), child_ids)
  chart_Ψ_bcast = lazy_map(Reindex(map(Broadcasting, Ψ_maps)), cell_to_chart)
  local_lazy    = lazy_map(evaluate, chart_Ψ_bcast, ref_per_fine)

  # ── 5. Resolve orientation_style ─────────────────────────────────────────
  os = if isnothing(orientation_style)
    Gridap.Geometry.OrientationStyle(Gridap.Geometry.get_grid(fine_model))
  else
    orientation_style
  end

  # ── 6. Build per-cell physical maps and metric lazily from cell_to_chart ─
  cell_phys_maps = lazy_map(Reindex(physical_maps), cell_to_chart)
  cell_metric    = lazy_map(Reindex(metric_fields),  cell_to_chart)

  atlas_grid = AtlasGrid(
    Gridap.Geometry.get_grid(fine_model),
    local_lazy,
    cell_phys_maps,
    cell_metric,
    os,
    manifold_style,
  )

  atlas_grid, fine_model
end

"""
    AtlasGrid(coarse_model, coarse_local_coords, physical_maps, num_refinements=1;
              metric_fields=nothing, orientation_style=nothing, manifold_style=ExtrinsicManifold())

Build an `AtlasGrid` from a coarse model.  When `metric_fields` is not provided it is
computed lazily as `Operation(JtJ())(gradient(φ))` for each physical map `φ`.
"""
function AtlasGrid(
    coarse_model        :: Gridap.Geometry.DiscreteModel{Dc,Dc},
    coarse_local_coords,
    physical_maps,
    num_refinements :: Int = 1;
    metric_fields     = nothing,
    orientation_style = nothing,
    manifold_style    = ExtrinsicManifold(),
) where Dc
  _metric = isnothing(metric_fields) ?
    [Operation(JtJ())(gradient(φ)) for φ in physical_maps] : metric_fields
  atlas_grid, _ = _build_atlas_grid(
    coarse_model, coarse_local_coords, physical_maps, _metric,
    num_refinements, orientation_style, manifold_style)
  atlas_grid
end

# ----------------------------------------------------------
# Gridap.Geometry.Grid{Dc,Dc} interface
# ----------------------------------------------------------

Gridap.Geometry.OrientationStyle(
  ::Type{<:AtlasGrid{Dc,Da,G,A,P,C,O,M}}) where {Dc,Da,G,A,P,C,O,M} = O()

Gridap.Geometry.num_cells(g::AtlasGrid) = length(g.cell_physical_maps)

Gridap.Geometry.get_reffes(g::AtlasGrid) = Gridap.Geometry.get_reffes(g.param_grid)

Gridap.Geometry.get_cell_type(g::AtlasGrid) = Gridap.Geometry.get_cell_type(g.param_grid)

Gridap.Geometry.get_cell_coordinates(g::AtlasGrid) = g.cell_local_coords

# Delegate to param_grid: gives the correct shared-node count for FESpace/num_nodes.
# Coordinate values are junk (2D parametric), but FEM assembly uses get_cell_map (overridden
# below) and never reads these values — only the length matters.
Gridap.Geometry.get_node_coordinates(g::AtlasGrid) =
  Gridap.Geometry.get_node_coordinates(g.param_grid)

# get_cell_map dispatches on ManifoldStyle:
#   ExtrinsicManifold — RefFE→Chart composed with Chart→Physical (current behavior)
#   IntrinsicManifold — RefFE→Chart only; metric queried via get_cell_metric
Gridap.Geometry.get_cell_map(g::AtlasGrid) = _get_cell_map(g, ManifoldStyle(g))

function _get_cell_map(g::AtlasGrid, ::ExtrinsicManifold)
  lazy_map(∘, g.cell_physical_maps, _chart_maps(g))
end

function _get_cell_map(g::AtlasGrid, ::IntrinsicManifold)
  _chart_maps(g)
end

# Delegate to param_grid: shared-node IDs consistent with get_node_coordinates.
# visualization_data builds its own DG sequential table for VTK.
Gridap.Geometry.get_cell_node_ids(g::AtlasGrid) =
  Gridap.Geometry.get_cell_node_ids(g.param_grid)

# ----------------------------------------------------------
# Custom API
# ----------------------------------------------------------

get_param_grid(g::AtlasGrid)                       = g.param_grid
get_cell_physical_maps(g::AtlasGrid)               = g.cell_physical_maps
get_cell_metric(g::AtlasGrid)                      = g.cell_metric
get_ambient_dim(::AtlasGrid{Dc,Da}) where {Dc,Da}  = Da

# ============================================================
# Coarse mesh library and convenience constructors
# ============================================================

include("CoarseMeshes.jl")

"""
    AtlasGrid(info::CoarseMeshInfo, num_refinements;
              orientation_style=nothing, manifold_style=ExtrinsicManifold())

Build an `AtlasGrid` from a `CoarseMeshInfo`, using the analytic metric fields
stored in `info.metric_fields`.
"""
function AtlasGrid(
    info            :: CoarseMeshInfo,
    num_refinements :: Int;
    orientation_style = nothing,
    manifold_style    = ExtrinsicManifold(),
)
  AtlasGrid(info.model, info.local_coords, info.physical_maps, num_refinements;
            metric_fields=info.metric_fields, orientation_style, manifold_style)
end

"""
    AtlasGrid(info::CoarseMeshInfo, custom_maps, num_refinements;
              orientation_style=nothing, manifold_style=ExtrinsicManifold())

Build an `AtlasGrid` from a `CoarseMeshInfo` with custom physical maps (overriding
`info.physical_maps`).  Metric fields are recomputed from `custom_maps` via `JtJ`.
"""
function AtlasGrid(
    info            :: CoarseMeshInfo,
    custom_maps,
    num_refinements :: Int;
    orientation_style = nothing,
    manifold_style    = ExtrinsicManifold(),
)
  custom_metrics = [Operation(JtJ())(gradient(φ)) for φ in custom_maps]
  AtlasGrid(info.model, info.local_coords, custom_maps, num_refinements;
            metric_fields=custom_metrics, orientation_style, manifold_style)
end

"""
    AtlasGrid(mesh::AbstractCoarseMesh, num_refinements;
              orientation_style=nothing, manifold_style=ExtrinsicManifold())

Build an `AtlasGrid` directly from a mesh descriptor (e.g. `CubedSphereMesh(1.0)`).
"""
function AtlasGrid(
    mesh            :: AbstractCoarseMesh,
    num_refinements :: Int;
    orientation_style = nothing,
    manifold_style    = ExtrinsicManifold(),
)
  AtlasGrid(get_coarse_mesh(mesh), num_refinements; orientation_style, manifold_style)
end
