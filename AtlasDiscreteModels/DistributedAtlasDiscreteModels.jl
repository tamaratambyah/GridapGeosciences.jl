# DistributedAtlasDiscreteModels.jl
#
# Distributed atlas discrete model built on p4est / GridapP4est.
# Extends AtlasGrid/AtlasDiscreteModel from AtlasGrids.jl/AtlasDiscreteModels.jl
# to the distributed setting via OctreeDistributedDiscreteModel.
#
# Key design change vs. old SBAtlasModels.jl:
#   - AtlasGrid now stores LOCAL (α,β) reference coords (not 3D physical coords).
#   - Physical Da-dimensional coords are computed only in visualization_data
#     (inherited from AtlasDiscreteModels.jl via _local_to_physical).
#   - The p4est infrastructure provides (α,β) coords per fine cell directly;
#     these are stored as-is into cell_local_coords.
#
# Run (e.g. 4 MPI processes):
#   mpiexec -n 4 julia --project=. <script>

include("AtlasDiscreteModels.jl")

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
# AtlasOctreeDistributedDiscreteModel
# ============================================================

"""
    AtlasOctreeDistributedDiscreteModel{A,B}

Distributed discrete model for an atlas-based 2D manifold mesh built on p4est.

# Fields
- `octree_dmodel` — underlying `OctreeDistributedDiscreteModel{2,2}` owning the
  p4est forest, MPI topology, and adaptive refinement support.
- `dmodel`        — `GenericDistributedDiscreteModel{2,2}` wrapping per-rank
  `AtlasDiscreteModel` instances, each carrying an `AtlasGrid{2,Da}` with local
  (α,β) reference coords and the per-chart physical maps.

The `DistributedDiscreteModel` interface delegates to `dmodel`.
Physical Da-dimensional coords are computed on demand in `visualization_data`
via `_local_to_physical` (defined in AtlasDiscreteModels.jl).
"""
struct AtlasOctreeDistributedDiscreteModel{
  A <: OctreeDistributedDiscreteModel{2,2},
  B <: GenericDistributedDiscreteModel{2,2},
} <: GridapDistributed.DistributedDiscreteModel{2,2}
  octree_dmodel :: A
  dmodel        :: B
end

# ----------------------------------------------------------
# DistributedDiscreteModel interface (delegate to dmodel)
# ----------------------------------------------------------

GridapDistributed.local_views(m::AtlasOctreeDistributedDiscreteModel) =
  local_views(m.dmodel)

GridapDistributed.get_cell_gids(m::AtlasOctreeDistributedDiscreteModel) =
  get_cell_gids(m.dmodel)

# ----------------------------------------------------------
# Custom API
# ----------------------------------------------------------

get_octree_dmodel(m::AtlasOctreeDistributedDiscreteModel) = m.octree_dmodel
get_dmodel(m::AtlasOctreeDistributedDiscreteModel)        = m.dmodel

# ============================================================
# Cubed-sphere constructor  (Dc = 2, Da = 3)
# ============================================================

"""
    AtlasOctreeDistributedDiscreteModel(ranks, radius; num_initial_uniform_refinements=0)

Build a distributed `AtlasOctreeDistributedDiscreteModel` for the standard cubed sphere:
6 QUAD panels with gnomonic projection from the 2D parametric square `[-π/4, π/4]²`
onto the sphere of the given `radius`.

The per-rank `AtlasGrid{2,3}` stores local (α,β) reference coords (not 3D physical
coords). Physical coords are computed lazily only during VTK visualization.

# Arguments
- `ranks`                           — MPI ranks array.
- `radius`                          — sphere radius (positive real).
- `num_initial_uniform_refinements` — uniform p4est refinements at construction
                                      (default 0 = 6-cell coarse mesh).
"""
function AtlasOctreeDistributedDiscreteModel(
    ranks,
    radius :: Real;
    num_initial_uniform_refinements :: Int = 0,
)
  # ── 6-panel QUAD coarse mesh in 2D parametric space ─────────────────────
  coarse_model = _create_parametric_octree_dmodel_coarse_model()

  # ── Coarse-cell corner (α,β) coordinates, one entry per panel ───────────
  coarse_cell_alpha_beta = setup_coarse_cell_vertices_alpha_beta_coordinates()

  # ── Panel ids for the 6 coarse cells (identity: panel i = tree i) ────────
  coarse_cell_panels = collect(1:NPANELS)

  # ── Physical maps: panel p → gnomonic projection (α,β) → (X,Y,Z) ────────
  physical_maps = [ForwardMap2D(p, Float64(radius)) for p in 1:NPANELS]

  # ── p4est refinement: get per-rank fine-cell (α,β) coords and panel ids ──
  octree_dmodel, cell_wise_alpha_beta_coords, cell_panels =
    _generate_octree_dmodel_alpha_beta_coordinates_and_panels(
      ranks,
      coarse_model,
      num_initial_uniform_refinements,
      coarse_cell_alpha_beta,
      coarse_cell_panels,
    )

  orientation_style = Gridap.Geometry.NonOriented()

  # ── Per-rank: build AtlasGrid with local (α,β) coords (no 3D conversion) ─
  atlas_models = map(
    local_views(octree_dmodel.dmodel),
    cell_wise_alpha_beta_coords,
    cell_panels,
  ) do omodel, alpha_beta_coords, panels_local

    param_grid    = Gridap.Geometry.get_grid(omodel)
    grid_topology = Gridap.Geometry.get_grid_topology(omodel)
    face_labeling = Gridap.Geometry.get_face_labeling(omodel)

    # alpha_beta_coords :: Table{Point{2,Float64}} — already the local reference
    # coords for each fine cell in chart panels_local[i].  Pass directly to
    # AtlasGrid; physical 3D coords are deferred to visualization_data.
    atlas_grid = AtlasGrid(
      param_grid,
      alpha_beta_coords,
      panels_local,
      physical_maps,
      orientation_style,
    )

    AtlasDiscreteModel(atlas_grid, grid_topology, face_labeling)
  end

  dmodel = GenericDistributedDiscreteModel(
    atlas_models, get_cell_gids(octree_dmodel.dmodel))

  AtlasOctreeDistributedDiscreteModel(octree_dmodel, dmodel)
end

# ============================================================
# Distributed test helper  (for future use)
# ============================================================

"""
    test_atlas_octree_model(model, radius)

Check that every local fine cell's physical corners lie on the unit sphere of the
given `radius`.  Call inside a `map(local_views(model.dmodel)) do local_model ...`
block or directly in a serial-compatible context.
"""
function test_atlas_octree_model(local_model :: AtlasDiscreteModel{Dc,Da}, radius::Real) where {Dc,Da}
  g     = local_model.atlas_grid
  phys  = _local_to_physical(g.cell_local_coords, g.cell_to_chart, g.physical_maps)
  ncells = length(phys)
  for i in 1:ncells
    for pt in phys[i]
      r = sqrt(sum(pt[k]^2 for k in 1:Da))
      @assert isapprox(r, radius; atol=1e-10) "cell $i: r=$r ≠ radius=$radius  pt=$pt"
    end
  end
  ncells
end
