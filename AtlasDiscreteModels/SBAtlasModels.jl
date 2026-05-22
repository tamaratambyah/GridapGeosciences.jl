using FillArrays
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Helpers

using GridapGeosciences
using GridapDistributed
using GridapDistributed: GenericDistributedDiscreteModel, get_cell_gids
using GridapP4est
using GridapP4est: OctreeDistributedDiscreteModel
using PartitionedArrays: local_views

import GridapGeosciences.Geometry: NPANELS
import GridapGeosciences.Fields: ForwardMap2D
import GridapGeosciences.Distributed: _create_parametric_octree_dmodel_coarse_model
import GridapGeosciences.Distributed: setup_coarse_cell_vertices_alpha_beta_coordinates
import GridapGeosciences.Distributed: _generate_octree_dmodel_alpha_beta_coordinates_and_panels

# ============================================================
# AtlasGrid
# ============================================================

"""
    AtlasGrid{Dc,Da,Tp,G,A,M,O} <: Gridap.Geometry.Grid{Dc,Da}

A DG-style grid representing a `Dc`-dimensional manifold embedded in a `Da`-dimensional
ambient space (e.g. Dc=2, Da=3 for a 2D sphere surface in ℝ³).

# Fields
- `topo_grid::G`         — Purely topological grid in `Dc`-dimensional parameter space
                           (coordinates are junk; only connectivity is used).
- `cell_coordinates::A`  — Per-cell ambient coordinates (Da-dimensional). DG-style:
                           each cell has its own vector of `Da`-dim corner points.
                           No shared nodes between cells.
- `cell_to_chart::Vector{Int}` — For each fine cell, the index of its coarse chart panel.
- `physical_maps::M`     — One chart map (Dc → Da) per coarse panel, used to produce
                           `cell_coordinates` from the parametric quadrant positions.
- `reffe`                — The single `LagrangianRefFE{Dc}` shared by all cells.
- `orientation_style::O` — `Gridap.Geometry.Oriented()` or `NonOriented()`.

# Notes
- `get_node_coordinates` is **not implemented**: there are no globally shared nodes.
  Use `get_cell_coordinates` to access the Da-dimensional corner positions per cell.
- `get_cell_map` is automatically derived from `get_cell_coordinates` (via the
  Gridap default implementation), returning Da-dimensional cell maps suitable for
  FEM integration on the manifold.
"""
struct AtlasGrid{Dc, Da, Tp,
                 G <: Gridap.Geometry.Grid,
                 A <: AbstractVector,
                 M,
                 O <: Gridap.Geometry.OrientationStyle} <: Gridap.Geometry.Grid{Dc,Da}
  topo_grid       :: G
  cell_coordinates :: A
  cell_to_chart   :: Vector{Int}
  physical_maps   :: M
  reffe           :: Gridap.ReferenceFEs.LagrangianRefFE{Dc}
  orientation_style :: O

  function AtlasGrid(
    topo_grid        :: Gridap.Geometry.Grid{Dc,Dc},
    cell_coordinates :: AbstractVector{<:AbstractVector{<:Point{Da,Tp}}},
    cell_to_chart    :: AbstractVector{<:Integer},
    physical_maps,
    reffe            :: Gridap.ReferenceFEs.LagrangianRefFE{Dc},
    orientation_style :: Gridap.Geometry.OrientationStyle
  ) where {Dc,Da,Tp}
    n = Gridap.Geometry.num_cells(topo_grid)
    @check length(cell_coordinates) == n
        "cell_coordinates length ($(length(cell_coordinates))) ≠ num_cells ($(n))"
    @check length(cell_to_chart) == n
        "cell_to_chart length ($(length(cell_to_chart))) ≠ num_cells ($(n))"
    G = typeof(topo_grid)
    A = typeof(cell_coordinates)
    M = typeof(physical_maps)
    O = typeof(orientation_style)
    new{Dc,Da,Tp,G,A,M,O}(
      topo_grid,
      cell_coordinates,
      Vector{Int}(cell_to_chart),
      physical_maps,
      reffe,
      orientation_style,
    )
  end
end

# ----------------------------------------------------------
# Gridap.Geometry.Grid{Dc,Da} interface
# ----------------------------------------------------------

# Encode orientation in the type so that the trait dispatches correctly.
Gridap.Geometry.OrientationStyle(
  ::Type{<:AtlasGrid{Dc,Da,Tp,G,A,M,O}}) where {Dc,Da,Tp,G,A,M,O} = O()

# Number of cells: derived from cell_to_chart (not from topo_grid.num_cells,
# so that ghost cells can be included if needed).
Gridap.Geometry.num_cells(g::AtlasGrid) = length(g.cell_to_chart)

# Single cell type (all cells are the same reffe).
Gridap.Geometry.get_reffes(g::AtlasGrid) =
  Fill(g.reffe, Gridap.Geometry.num_cells(g))

Gridap.Geometry.get_cell_type(g::AtlasGrid) =
  Fill(Int8(1), Gridap.Geometry.num_cells(g))

# Connectivity comes from the underlying topological grid.
Gridap.Geometry.get_cell_node_ids(g::AtlasGrid) =
  Gridap.Geometry.get_cell_node_ids(g.topo_grid)

# DG-style: no shared nodes. Raise an error if accidentally called.
function Gridap.Geometry.get_node_coordinates(g::AtlasGrid)
  @notimplemented """\n
  AtlasGrid is DG-style: there are no globally shared nodes.
  Use `get_cell_coordinates(g)` to access the $(num_point_dims(g))-dimensional
  corner coordinates per cell instead.
  """
end

# Override to return the pre-computed Da-dimensional cell coordinates directly.
# This also makes `get_cell_map` (which calls `get_cell_coordinates` by default)
# return Da-dimensional maps — no further override needed.
Gridap.Geometry.get_cell_coordinates(g::AtlasGrid) = g.cell_coordinates

# ----------------------------------------------------------
# Custom API (not part of the standard Gridap interface)
# ----------------------------------------------------------

"""Return the underlying topological grid (Dc-dimensional, junk coordinates)."""
get_topo_grid(g::AtlasGrid) = g.topo_grid

"""Return the chart maps (one per coarse panel): Dc → Da."""
get_physical_maps(g::AtlasGrid) = g.physical_maps

"""Return the fine-cell → coarse-chart-id mapping."""
get_cell_to_chart(g::AtlasGrid) = g.cell_to_chart

"""Return the ambient (embedding) space dimension Da."""
get_ambient_dim(::AtlasGrid{Dc,Da}) where {Dc,Da} = Da


# ============================================================
# AtlasDiscreteModel
# ============================================================

"""
    AtlasDiscreteModel{Dc,Da,Tp,G,A,M,O,T,L} <: Gridap.Geometry.DiscreteModel{Dc,Dc}

A discrete model combining:
- an `AtlasGrid{Dc,Da}` for Da-dimensional manifold geometry (DG-style),
- a `GridTopology{Dc,Dc}` for the Dc-dimensional parametric topology (from p4est),
- a `FaceLabeling` for boundary/interface tags.

The model is typed as `DiscreteModel{Dc,Dc}` so that the standard Gridap
topology machinery (adjacency, face labeling) operates in parametric space.
The Da-dimensional geometry is accessed through the `atlas_grid` field.

# Notes
- `get_grid(model)` returns the `AtlasGrid{Dc,Da}`, not a standard
  `UnstructuredGrid`. Callers needing Da-dimensional cell maps should use
  `get_cell_map(get_grid(model))` or `get_cell_coordinates(get_grid(model))`.
- `get_grid_topology(model)` returns the Dc-dimensional topology (from p4est).
"""
struct AtlasDiscreteModel{Dc, Da, Tp, G, A, M, O,
                           T <: Gridap.Geometry.GridTopology,
                           L <: Gridap.Geometry.FaceLabeling
                           } <: Gridap.Geometry.DiscreteModel{Dc,Dc}
  atlas_grid    :: AtlasGrid{Dc,Da,Tp,G,A,M,O}
  grid_topology :: T
  face_labeling :: L

  function AtlasDiscreteModel(
    atlas_grid    :: AtlasGrid{Dc,Da,Tp,G,A,M,O},
    grid_topology :: Gridap.Geometry.GridTopology{Dc,Dc},
    face_labeling :: Gridap.Geometry.FaceLabeling,
  ) where {Dc,Da,Tp,G,A,M,O}
    T = typeof(grid_topology)
    L = typeof(face_labeling)
    new{Dc,Da,Tp,G,A,M,O,T,L}(atlas_grid, grid_topology, face_labeling)
  end
end

# ----------------------------------------------------------
# Gridap.Geometry.DiscreteModel{Dc,Dc} interface
# ----------------------------------------------------------

Gridap.Geometry.get_grid(m::AtlasDiscreteModel)          = m.atlas_grid
Gridap.Geometry.get_grid_topology(m::AtlasDiscreteModel) = m.grid_topology
Gridap.Geometry.get_face_labeling(m::AtlasDiscreteModel) = m.face_labeling

# ----------------------------------------------------------
# Custom API
# ----------------------------------------------------------

"""Return the AtlasGrid{Dc,Da} (same as get_grid, for clarity)."""
get_atlas_grid(m::AtlasDiscreteModel) = m.atlas_grid

"""Return the ambient (embedding) space dimension Da."""
get_ambient_dim(m::AtlasDiscreteModel{Dc,Da}) where {Dc,Da} = Da

"""Return the chart maps (one per coarse panel): Dc → Da."""
get_physical_maps(m::AtlasDiscreteModel) = get_physical_maps(m.atlas_grid)

"""Return the fine-cell → coarse-chart-id mapping."""
get_cell_to_chart(m::AtlasDiscreteModel) = get_cell_to_chart(m.atlas_grid)


# ============================================================
# Helper: transform 2D (α,β) → 3D ambient coordinates
# ============================================================

"""
    _transform_alpha_beta_to_3d(alpha_beta_table, panels_local, physical_maps)

Transform a `Table` of 2D parametric `(α,β)` cell corner coordinates to `Da`-dimensional
ambient coordinates by applying the per-panel chart map (forward map Dc → Da).

# Arguments
- `alpha_beta_table` : `Table{Point{2,Float64}}` — one sub-array of 2D corners per cell.
- `panels_local`     : `AbstractVector{Int}` — chart/panel index (1-based) for each cell.
- `physical_maps`    : indexable collection of `ForwardMap2D`, indexed by panel id.

# Returns
`Gridap.Arrays.Table{Point{3,Float64}}` with the same pointer structure as
`alpha_beta_table` but with 3D corner coordinates.
"""
function _transform_alpha_beta_to_3d(
    alpha_beta_table,
    panels_local::AbstractVector{Int},
    physical_maps,
)
  n_data = length(alpha_beta_table.data)
  data_3d = Vector{Point{3,Float64}}(undef, n_data)
  for i in 1:length(panels_local)
    fwd_map       = physical_maps[panels_local[i]]
    start         = alpha_beta_table.ptrs[i]
    stop          = alpha_beta_table.ptrs[i+1] - 1
    coords_2d     = alpha_beta_table.data[start:stop]   # Vector{Point{2,Float64}}, copy
    cache         = Gridap.Arrays.return_cache(fwd_map, coords_2d)
    coords_3d     = Gridap.Arrays.evaluate!(cache, fwd_map, coords_2d)
    data_3d[start:stop] = coords_3d
  end
  Gridap.Arrays.Table(data_3d, copy(alpha_beta_table.ptrs))
end

# ============================================================
# AtlasOctreeDistributedDiscreteModel
# ============================================================

"""
    AtlasOctreeDistributedDiscreteModel{A,B} <: GridapDistributed.DistributedDiscreteModel{2,2}

Distributed discrete model for an atlas-based 2D manifold mesh built on p4est.

# Fields
- `octree_dmodel` : the underlying `OctreeDistributedDiscreteModel{2,2}` that owns the
  p4est forest, the MPI topology, and supports adaptive refinement.
- `dmodel`        : a `GenericDistributedDiscreteModel{2,2}` wrapping per-rank
  `AtlasDiscreteModel` instances (each carries an `AtlasGrid{2,Da}` with Da-dimensional
  cell coordinates and the chart maps).

The `DistributedDiscreteModel` interface is delegated to `dmodel`.
"""
struct AtlasOctreeDistributedDiscreteModel{
  A <: OctreeDistributedDiscreteModel{2,2},
  B <: GenericDistributedDiscreteModel{2,2},
} <: GridapDistributed.DistributedDiscreteModel{2,2}
  octree_dmodel :: A
  dmodel        :: B
end

# ----------------------------------------------------------
# GridapDistributed.DistributedDiscreteModel interface (delegate to dmodel)
# ----------------------------------------------------------

GridapDistributed.local_views(m::AtlasOctreeDistributedDiscreteModel) =
  local_views(m.dmodel)

GridapDistributed.get_cell_gids(m::AtlasOctreeDistributedDiscreteModel) =
  get_cell_gids(m.dmodel)

# ----------------------------------------------------------
# Custom API
# ----------------------------------------------------------

"""Return the underlying OctreeDistributedDiscreteModel (p4est owner)."""
get_octree_dmodel(m::AtlasOctreeDistributedDiscreteModel) = m.octree_dmodel

"""Return the GenericDistributedDiscreteModel of per-rank AtlasDiscreteModels."""
get_dmodel(m::AtlasOctreeDistributedDiscreteModel) = m.dmodel

# ============================================================
# Cubed-sphere constructor  (Dc = 2, Da = 3)
# ============================================================

"""
    AtlasOctreeDistributedDiscreteModel(ranks, radius; num_initial_uniform_refinements=0)

Build a distributed `AtlasOctreeDistributedDiscreteModel` for the standard cubed sphere:
6 QUAD panels with gnomonic projection from the 2D parametric square `[-π/4, π/4]²`
onto the sphere of the given `radius`.

# Arguments
- `ranks`                           : MPI ranks array (e.g. `distribute_with_mpi(LinearIndices((np,)))`).
- `radius`                          : sphere radius (positive real number).
- `num_initial_uniform_refinements` : number of uniform p4est refinements applied at construction
                                      (default `0` gives the 6-cell coarse mesh).

# Returns
An `AtlasOctreeDistributedDiscreteModel` whose per-rank `AtlasDiscreteModel` instances
each hold an `AtlasGrid{2,3}` with 3D (sphere surface) cell corner coordinates.
"""
function AtlasOctreeDistributedDiscreteModel(
    ranks,
    radius::Real;
    num_initial_uniform_refinements::Int = 0,
)
  # ── 6-panel QUAD coarse mesh in 2D parametric space (junk node coords) ──────
  coarse_model = _create_parametric_octree_dmodel_coarse_model()

  # ── Coarse-cell corner coordinates in (α,β) ∈ [-π/4,π/4]² per panel ────────
  coarse_cell_alpha_beta = setup_coarse_cell_vertices_alpha_beta_coordinates()

  # ── One-based panel ids for the 6 coarse cells (identity: panel i → tree i) ─
  coarse_cell_panels = collect(1:NPANELS)

  # ── Physical maps: panel p → gnomonic projection Dc=2 → Da=3 ────────────────
  physical_maps = [ForwardMap2D(p, Float64(radius)) for p in 1:NPANELS]

  # ── Build OctreeDistributedDiscreteModel + per-rank 2D cell corners + panels ─
  octree_dmodel, cell_wise_alpha_beta_coords, cell_panels =
    _generate_octree_dmodel_alpha_beta_coordinates_and_panels(
      ranks,
      coarse_model,
      num_initial_uniform_refinements,
      coarse_cell_alpha_beta,
      coarse_cell_panels,
    )

  reffe = Gridap.ReferenceFEs.ReferenceFE(
    Gridap.ReferenceFEs.QUAD, Gridap.ReferenceFEs.lagrangian, Float64, 1)
  orientation_style = Gridap.Geometry.NonOriented()

  # ── Per-rank: transform 2D (α,β) → 3D sphere coords; build AtlasDiscreteModel
  atlas_models = map(
    local_views(octree_dmodel.dmodel),
    cell_wise_alpha_beta_coords,
    cell_panels,
  ) do omodel, alpha_beta_coords, panels_local

    cell_coords_3d = _transform_alpha_beta_to_3d(
      alpha_beta_coords, panels_local, physical_maps)

    topo_grid     = Gridap.Geometry.get_grid(omodel)
    grid_topology = Gridap.Geometry.get_grid_topology(omodel)
    face_labeling = Gridap.Geometry.get_face_labeling(omodel)

    atlas_grid = AtlasGrid(
      topo_grid,
      cell_coords_3d,
      panels_local,
      physical_maps,
      reffe,
      orientation_style,
    )

    AtlasDiscreteModel(atlas_grid, grid_topology, face_labeling)
  end

  dmodel = GenericDistributedDiscreteModel(
    atlas_models, get_cell_gids(octree_dmodel.dmodel))

  AtlasOctreeDistributedDiscreteModel(octree_dmodel, dmodel)
end
