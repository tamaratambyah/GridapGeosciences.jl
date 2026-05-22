# test_SBAtlasModels.jl
#
# Verify that AtlasOctreeDistributedDiscreteModel (SBAtlasModels.jl) produces
# correct 3D sphere-surface cell coordinates.
#
# Strategy: build both models with the same refinement level, then for every
# fine cell i check:
#
#   actual_3d[i]   = cell_coordinates(new AtlasGrid{2,3})[i]
#   expected_3d[i] = ForwardMap2D(panel[i], radius).(old_2d_coords[i])
#
# where old_2d_coords come from the cell_map of the reference
# CubedSphere2DParametricOctreeDistributedDiscreteModel, which stores exactly
# the (α,β) Table that we also feed into _transform_alpha_beta_to_3d.
#
# Run (single process):
#   julia --project=. test_SBAtlasModels.jl
# or with MPI (e.g. 4 processes):
#   mpiexec -n 4 julia --project=. test_SBAtlasModels.jl

using FillArrays
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Helpers
using GridapGeosciences
using GridapDistributed
using GridapP4est
using PartitionedArrays
using MPI

import GridapGeosciences.Fields: ForwardMap2D

include("SBAtlasModels.jl")

# ── MPI setup ─────────────────────────────────────────────────────────────────
MPI.Init()
nprocs = MPI.Comm_size(MPI.COMM_WORLD)
ranks  = distribute_with_mpi(LinearIndices((nprocs,)))

const RADIUS  = 1.0
const NUM_REF = 1     # uniform refinements (coarse mesh = 6 cells, level 1 = 24 cells)

# ── Build models ──────────────────────────────────────────────────────────────
# New model: AtlasGrid{2,3} with 3D sphere-surface cell coordinates
model_new = AtlasOctreeDistributedDiscreteModel(
  ranks, RADIUS; num_initial_uniform_refinements = NUM_REF)

# Reference model: CubedSphere2DParametricDiscreteModel with 2D (α,β) cell coords
model_old = CubedSphere2DParametricOctreeDistributedDiscreteModel(
  ranks, RADIUS; num_initial_uniform_refinements = NUM_REF)

# ── Per-rank comparison ───────────────────────────────────────────────────────
map(
  local_views(model_new.dmodel),
  local_views(model_old.parametric_dmodel),
) do local_new, local_old

  # --- new model side ---
  # cell_coordinates :: Table{Point{3,Float64}}  (3D sphere surface coords)
  # cell_to_chart    :: Vector{Int}              (panel id per cell, 1-based)
  cell_coords_3d = local_new.atlas_grid.cell_coordinates
  panels         = local_new.atlas_grid.cell_to_chart

  # --- old model side ---
  # local_old :: CubedSphere2DParametricDiscreteModel
  # local_old.grid.cell_map = lazy_map(linear_combination, alpha_beta_table, shape_funs)
  # .args[1] is the Table{Point{2,Float64}} of (α,β) corners fed to linear_combination.
  cell_coords_2d = local_old.grid.cell_map.args[1]
  panel_ids_old  = local_old.panel_ids

  @assert length(cell_coords_3d) == length(cell_coords_2d) string(
    "Cell count mismatch: new=$(length(cell_coords_3d)), old=$(length(cell_coords_2d))")

  @assert panels == Vector{Int}(panel_ids_old) string(
    "Panel id arrays differ between new and old model")

  ncells = length(panels)
  for i in 1:ncells
    fwd_map    = ForwardMap2D(panels[i], RADIUS)
    coords_2d  = cell_coords_2d[i]                  # Vector{Point{2,Float64}}
    cache      = Gridap.Arrays.return_cache(fwd_map, coords_2d)
    expected   = Gridap.Arrays.evaluate!(cache, fwd_map, coords_2d)   # Vector{Point{3}}
    actual     = cell_coords_3d[i]                  # Vector{Point{3,Float64}}

    for (k, (e, a)) in enumerate(zip(expected, actual))
      @assert e ≈ a string(
        "Rank $(MPI.Comm_rank(MPI.COMM_WORLD)) cell $i corner $k ",
        "panel=$(panels[i]): expected=$e  got=$a")
    end
  end

  rank = MPI.Comm_rank(MPI.COMM_WORLD)
  println("Rank $rank: $ncells cells verified ✓")
end

println("test_SBAtlasModels: ALL CHECKS PASSED")
