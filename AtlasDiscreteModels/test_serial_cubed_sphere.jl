# test_serial_cubed_sphere.jl
#
# Tests AtlasGrid / AtlasDiscreteModel on the cubed sphere (6 QUAD panels,
# gnomonic projection).  Serial — no MPI required.
#
# Checks:
#   1. Cell count = NPANELS × 4^NUM_REF
#   2. Each panel has exactly 4^NUM_REF fine cells (cell_to_chart)
#   3. Local (α,β) coords are 2D and lie in [-π/4, π/4]²
#   4. Physical coords lie on the sphere: ‖(X,Y,Z)‖ ≈ RADIUS
#   5. writevtk produces a 3D sphere surface VTK file
#
# Run:
#   julia --project=. AtlasDiscreteModels/test_serial_cubed_sphere.jl

using Gridap
using Gridap.Geometry, Gridap.Arrays, Gridap.ReferenceFEs
using GridapGeosciences

include("AtlasDiscreteModels.jl")

import GridapGeosciences.Geometry: NPANELS
import GridapGeosciences.Fields: ForwardMap2D
import GridapGeosciences.Distributed: _create_parametric_octree_dmodel_coarse_model
import GridapGeosciences.Distributed: setup_coarse_cell_vertices_alpha_beta_coordinates

const RADIUS  = 1.0
const NUM_REF = 2    # 6 coarse cells → 6 × 4² = 96 fine cells

# ── Coarse 6-panel QUAD mesh (correct topology, junk coordinates) ─────────────
coarse_model = _create_parametric_octree_dmodel_coarse_model()

# ── Coarse-cell corner (α,β) ∈ [-π/4, π/4]² — same corners for all 6 panels ──
coarse_local_coords = setup_coarse_cell_vertices_alpha_beta_coordinates()

# ── Physical maps: gnomonic projection panel p → sphere surface ───────────────
physical_maps = [ForwardMap2D(p, Float64(RADIUS)) for p in 1:NPANELS]

# ── Build AtlasDiscreteModel ──────────────────────────────────────────────────
atlas_model = AtlasDiscreteModel(coarse_model, coarse_local_coords, physical_maps, NUM_REF)
atlas_grid  = get_atlas_grid(atlas_model)

println("AtlasGrid:  Dc=$(Gridap.Geometry.num_cell_dims(atlas_grid)), Da=$(get_ambient_dim(atlas_grid))")
println("            cells = $(Gridap.Geometry.num_cells(atlas_grid))")
println("AtlasModel: built ✓")

# ── 1. Cell count ─────────────────────────────────────────────────────────────
expected_cells = NPANELS * 4^NUM_REF
@assert Gridap.Geometry.num_cells(atlas_grid) == expected_cells "Expected $expected_cells cells, got $(Gridap.Geometry.num_cells(atlas_grid))"

# ── 2. All 6 panels present, each with exactly 4^NUM_REF fine cells ───────────
c2c = get_cell_to_chart(atlas_grid)
cells_per_panel = 4^NUM_REF
for p in 1:NPANELS
  count = sum(c2c .== p)
  @assert count == cells_per_panel "Panel $p: expected $cells_per_panel cells, got $count"
end
println("Panel assignment check passed ✓")

# ── 3. Local (α,β) coords are 2D and lie in [-π/4, π/4]² ─────────────────────
local_coords = Gridap.Geometry.get_cell_coordinates(atlas_grid)
half_edge = π / 4
for i in 1:Gridap.Geometry.num_cells(atlas_grid)
  for pt in local_coords[i]
    @assert -half_edge - 1e-12 ≤ pt[1] ≤ half_edge + 1e-12 "cell $i: α=$(pt[1]) out of [-π/4, π/4]"
    @assert -half_edge - 1e-12 ≤ pt[2] ≤ half_edge + 1e-12 "cell $i: β=$(pt[2]) out of [-π/4, π/4]"
  end
end
println("Local coord range check passed ✓")

# ── 4. Physical coords lie on the sphere ──────────────────────────────────────
phys_coords = _local_to_physical(
  atlas_grid.cell_local_coords,
  atlas_grid.cell_to_chart,
  atlas_grid.physical_maps,
)
for i in 1:Gridap.Geometry.num_cells(atlas_grid)
  for pt in phys_coords[i]
    r = sqrt(pt[1]^2 + pt[2]^2 + pt[3]^2)
    @assert isapprox(r, RADIUS; atol=1e-12) "cell $i: r=$r ≠ RADIUS=$RADIUS  pt=$pt"
  end
end
println("Physical coord sphere check passed ✓")

# ── 5. VTK output ─────────────────────────────────────────────────────────────
mkpath("output")
writevtk(atlas_model, "output/cubed_sphere_serial")
println("Written output/cubed_sphere_serial_2.vtu — open in Paraview ✓")

println("test_serial_cubed_sphere: ALL CHECKS PASSED")
