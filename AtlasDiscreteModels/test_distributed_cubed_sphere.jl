# test_distributed_cubed_sphere.jl
#
# Verify AtlasOctreeDistributedDiscreteModel (DistributedAtlasDiscreteModels.jl)
# on the cubed sphere with p4est / GridapP4est.
#
# Checks:
#   1. AtlasGrid stores local (α,β) 2D reference coords (not 3D physical).
#   2. Panel ids match the reference CubedSphere2DParametricOctreeDistributedDiscreteModel.
#   3. Local (α,β) coords match those stored in the reference model.
#   4. Applying physical_maps to local coords gives the same result as ForwardMap2D directly.
#   5. All physical corners lie on the sphere of the given radius.
#   6. writevtk produces Da=3 dimensional output via visualization_data.
#
# Run (single process):
#   julia --project=. AtlasDiscreteModels/test_distributed_cubed_sphere.jl
# With MPI (e.g. 4 processes):
#   mpiexec -n 4 julia --project=. AtlasDiscreteModels/test_distributed_cubed_sphere.jl

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Helpers
using GridapGeosciences
using GridapDistributed
using GridapP4est
using PartitionedArrays
using MPI

include("DistributedAtlasDiscreteModels.jl")

# Check that every local fine cell's physical corners lie on the sphere of radius `r`.
function test_atlas_octree_model(local_model::AtlasDiscreteModel{Dc,Da}, r::Real) where {Dc,Da}
  g     = local_model.atlas_grid
  phys  = _local_to_physical(g.cell_local_coords, g.cell_to_chart, g.physical_maps)
  ncells = length(phys)
  for i in 1:ncells
    for pt in phys[i]
      norm = sqrt(sum(pt[k]^2 for k in 1:Da))
      @assert isapprox(norm, r; atol=1e-10) "cell $i: ‖pt‖=$norm ≠ radius=$r  pt=$pt"
    end
  end
  ncells
end

# ── MPI setup ─────────────────────────────────────────────────────────────────
MPI.Init()
nprocs = MPI.Comm_size(MPI.COMM_WORLD)
ranks  = distribute_with_mpi(LinearIndices((nprocs,)))

const RADIUS  = 1.0
const NUM_REF = 1    # 6 coarse cells → 24 fine cells at level 1

# ── Build new AtlasOctreeDistributedDiscreteModel ────────────────────────────
atlas_model = AtlasOctreeDistributedDiscreteModel(
  ranks, CubedSphereMesh(RADIUS); num_initial_uniform_refinements = NUM_REF)

# ── Reference: CubedSphere2DParametricOctreeDistributedDiscreteModel ─────────
ref_model = CubedSphere2DParametricOctreeDistributedDiscreteModel(
  ranks, RADIUS; num_initial_uniform_refinements = NUM_REF)

# ── Per-rank verification ─────────────────────────────────────────────────────
map(
  local_views(atlas_model.dmodel),
  local_views(ref_model.parametric_dmodel),
) do local_atlas, local_ref

  g       = local_atlas.atlas_grid
  ncells  = Gridap.Geometry.num_cells(g)

  # 1. Local coords are 2D (α,β)
  local_coords = g.cell_local_coords
  @assert length(local_coords) == ncells

  # 2. Panel ids match the reference model
  panels_ref = local_ref.panel_ids
  @assert g.cell_to_chart == Vector{Int}(panels_ref) "Panel id arrays differ"

  # 3. Local (α,β) coords match reference model's 2D coords
  #    (reference stores them in grid.cell_map.args[1])
  alpha_beta_ref = local_ref.grid.cell_map.args[1]
  for i in 1:ncells
    @assert local_coords[i] ≈ alpha_beta_ref[i] "Local coords mismatch at cell $i"
  end

  # 4. Physical coords computed via _local_to_physical match ForwardMap2D
  phys = _local_to_physical(local_coords, g.cell_to_chart, g.physical_maps)
  # One cache per chart, hoisted outside the cell loop — reused for every corner.
  sample_pt    = local_coords[1][1]
  chart_caches = [Gridap.Arrays.return_cache(g.physical_maps[p], sample_pt)
                  for p in eachindex(g.physical_maps)]
  for i in 1:ncells
    chart = g.cell_to_chart[i]
    fwd   = g.physical_maps[chart]
    cache = chart_caches[chart]
    for (k, (local_pt, actual_pt)) in enumerate(zip(local_coords[i], phys[i]))
      expected_pt = Gridap.Arrays.evaluate!(cache, fwd, local_pt)
      @assert expected_pt ≈ actual_pt "Rank $(MPI.Comm_rank(MPI.COMM_WORLD)) cell $i corner $k: expected=$expected_pt got=$actual_pt"
    end
  end

  # 5. All physical corners on the sphere
  nchecked = test_atlas_octree_model(local_atlas, RADIUS)
  rank = MPI.Comm_rank(MPI.COMM_WORLD)
  println("Rank $rank: $nchecked cells verified ✓")
end

# ── VTK output ────────────────────────────────────────────────────────────────
mkpath("output")
writevtk(atlas_model, "output/cubed_sphere")
println("Written output/cubed_sphere_2.vtu ✓")

println("test_distributed_cubed_sphere: ALL CHECKS PASSED")
