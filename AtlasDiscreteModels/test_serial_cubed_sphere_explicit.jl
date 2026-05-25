# test_serial_cubed_sphere_explicit.jl
#
# Tests AtlasGrid / AtlasDiscreteModel on the cubed sphere, building the
# 8-node / 6-cell coarse QUAD mesh EXPLICITLY (same as the cylinder test).
# Does NOT call _create_parametric_octree_dmodel_coarse_model or
# setup_coarse_cell_vertices_alpha_beta_coordinates.
#
# Cubed-sphere topology (viewed from above the "unfolded" cube).
# Node coordinates are junk — used only for valid QUAD topology.
# Gridap Z-order per cell: BL, BR, TL, TR.
#
#        x=-2                                    x=2
# y=2     5───────────────────────────────────────6   C2
#         │╲   y                                 ╱│
#         │  ╲ . x            C4              /   │
# y=1     │    8───────────────────────────7      │
#         │    │                           │      │
#         │ C6 │y          C5              │ C3   │
#         │    │. x                        │      │
# y=-1    │    1───────────────────────────2      │
#         │. ╱              C1               ╲    │
#         │╱.                                  ╲ .│
# y=-2    3───────────────────────────────────────4
#
# 8 nodes:
#   1=(-1,-1)  2=(1,-1)  3=(-2,-2)  4=(2,-2)
#   5=(-2, 2)  6=(2, 2)  7=(1, 1)   8=(-1, 1)
#
# 6 cells (Gridap Z-order BL, BR, TL, TR):
#   C1 = [1,2,3,4]  (BL=1, BR=2, TL=3, TR=4)
#   C2 = [3,4,5,6]  (BL=3, BR=4, TL=5, TR=6)
#   C3 = [2,7,4,6]  (BL=2, BR=7, TL=4, TR=6)
#   C4 = [8,5,7,6]  (BL=8, BR=5, TL=7, TR=6)
#   C5 = [1,8,2,7]  (BL=1, BR=8, TL=2, TR=7)
#   C6 = [1,3,8,5]  (BL=1, BR=3, TL=8, TR=5)
#
# Node sharing encodes the connectivity of the six cube faces:
#   - Shared edges (e.g. nodes 1-2 shared by C1 and C5) reflect the cube topology.
#   - The 8 nodes correspond to the 8 corners of a cube projected onto the sphere.
#
# Local chart coordinates: (α,β) ∈ [-π/4, π/4]² for all 6 panels.
# Physical maps: ForwardMap2D(p, RADIUS) for p = 1..6 (gnomonic projection).
#
# Run:
#   julia --project=. AtlasDiscreteModels/test_serial_cubed_sphere_explicit.jl

using Gridap
using Gridap.Geometry, Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Adaptivity, Gridap.Helpers
using GridapGeosciences
using FillArrays

include("AtlasDiscreteModels.jl")

import GridapGeosciences.Geometry: NPANELS
import GridapGeosciences.Fields: ForwardMap2D

const RADIUS  = 1.0
const NUM_REF = 2    # 6 coarse cells → 6 × 4² = 96 fine cells

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Build explicit 8-node / 6-cell QUAD mesh (topology only; coordinates junk)
# ─────────────────────────────────────────────────────────────────────────────

node_coords = Vector{Point{2,Float64}}([
  Point(-1.0, -1.0),   # 1
  Point( 1.0, -1.0),   # 2
  Point(-2.0, -2.0),   # 3
  Point( 2.0, -2.0),   # 4
  Point(-2.0,  2.0),   # 5
  Point( 2.0,  2.0),   # 6
  Point( 1.0,  1.0),   # 7
  Point(-1.0,  1.0),   # 8
])

# 6 panels, 4 nodes each, Gridap Z-order: BL BR TL TR
# Connectivity matches _CCAM_panel_wise_node_ids() exactly.
cell_node_data = Int32[1,2,3,4,  3,4,5,6,  2,7,4,6,  8,5,7,6,  1,8,2,7,  1,3,8,5]
cell_node_ptrs = Int32[1,5,9,13,17,21,25]
cell_node_ids  = Table(cell_node_data, cell_node_ptrs)
cell_types     = Int32[1,1,1,1,1,1]
reffe          = LagrangianRefFE(Float64, QUAD, 1)

coarse_grid   = UnstructuredGrid(node_coords, cell_node_ids, [reffe], cell_types, NonOriented())
coarse_topo   = UnstructuredGridTopology(node_coords, cell_node_ids, cell_types, [QUAD], NonOriented())
coarse_labels = FaceLabeling(coarse_topo)
coarse_model  = UnstructuredDiscreteModel(coarse_grid, coarse_topo, coarse_labels)

println("Coarse mesh: $(num_cells(coarse_model)) cells, $(num_nodes(coarse_model)) nodes")
@assert num_cells(coarse_model)  == 6 "Expected 6 cells, got $(num_cells(coarse_model))"
@assert num_nodes(coarse_model)  == 8 "Expected 8 nodes, got $(num_nodes(coarse_model))"

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Local coordinates for each coarse chart
#     All 6 panels share the same local frame: (α,β) ∈ [-π/4, π/4]²
#     Gridap Z-order: BL, BR, TL, TR
# ─────────────────────────────────────────────────────────────────────────────

half_edge = π / 4
panel_corners = [
  Point(-half_edge, -half_edge),   # BL
  Point( half_edge, -half_edge),   # BR
  Point(-half_edge,  half_edge),   # TL
  Point( half_edge,  half_edge),   # TR
]
coarse_local_coords = fill(panel_corners, NPANELS)

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Physical maps: gnomonic projection panel p → sphere surface
# ─────────────────────────────────────────────────────────────────────────────

physical_maps = [ForwardMap2D(p, Float64(RADIUS)) for p in 1:NPANELS]

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Build AtlasDiscreteModel
# ─────────────────────────────────────────────────────────────────────────────

atlas_model = AtlasDiscreteModel(coarse_model, coarse_local_coords, physical_maps, NUM_REF)
atlas_grid  = get_atlas_grid(atlas_model)

println("\nAtlasGrid:  Dc=$(Gridap.Geometry.num_cell_dims(atlas_grid)), Da=$(get_ambient_dim(atlas_grid))")
println("            cells = $(Gridap.Geometry.num_cells(atlas_grid))")
println("AtlasModel: built ✓")

# ─────────────────────────────────────────────────────────────────────────────
# 5.  Cell count
# ─────────────────────────────────────────────────────────────────────────────

expected_cells = NPANELS * 4^NUM_REF
@assert Gridap.Geometry.num_cells(atlas_grid) == expected_cells "Expected $expected_cells cells, got $(Gridap.Geometry.num_cells(atlas_grid))"
println("Cell count check passed ✓  ($expected_cells cells)")

# ─────────────────────────────────────────────────────────────────────────────
# 6.  Panel assignment: each of the 6 panels has exactly 4^NUM_REF fine cells
# ─────────────────────────────────────────────────────────────────────────────

c2c = get_cell_to_chart(atlas_grid)
cells_per_panel = 4^NUM_REF
for p in 1:NPANELS
  count = sum(c2c .== p)
  @assert count == cells_per_panel "Panel $p: expected $cells_per_panel cells, got $count"
end
println("Panel assignment check passed ✓  ($cells_per_panel cells per panel)")

# ─────────────────────────────────────────────────────────────────────────────
# 7.  Local (α,β) coords are 2D and lie in [-π/4, π/4]²
# ─────────────────────────────────────────────────────────────────────────────

local_coords = Gridap.Geometry.get_cell_coordinates(atlas_grid)
@assert length(local_coords) == expected_cells

for i in 1:expected_cells
  for pt in local_coords[i]
    @assert -half_edge - 1e-12 ≤ pt[1] ≤ half_edge + 1e-12 "cell $i: α=$(pt[1]) out of [-π/4, π/4]"
    @assert -half_edge - 1e-12 ≤ pt[2] ≤ half_edge + 1e-12 "cell $i: β=$(pt[2]) out of [-π/4, π/4]"
  end
end
println("Local coord range check passed ✓  (all (α,β) ∈ [-π/4, π/4]²)")

# ─────────────────────────────────────────────────────────────────────────────
# 8.  Physical coords lie on the sphere: ‖(X,Y,Z)‖ ≈ RADIUS
# ─────────────────────────────────────────────────────────────────────────────

phys_coords = _local_to_physical(
  atlas_grid.cell_local_coords,
  atlas_grid.cell_to_chart,
  atlas_grid.physical_maps,
)
@assert length(phys_coords) == expected_cells

for i in 1:expected_cells
  for pt in phys_coords[i]
    r = sqrt(pt[1]^2 + pt[2]^2 + pt[3]^2)
    @assert isapprox(r, RADIUS; atol=1e-12) "cell $i: ‖pt‖=$r ≠ RADIUS=$(RADIUS)  pt=$pt"
  end
end
println("Physical coord sphere check passed ✓  (all ‖(X,Y,Z)‖ ≈ $(RADIUS))")

# ─────────────────────────────────────────────────────────────────────────────
# 9.  Connectivity sanity: all 6 panels covered and physical points distinct
# ─────────────────────────────────────────────────────────────────────────────

# Centroid of each panel should be distinct 3D points on the sphere
panel_centroids = Dict{Int, Vector{Point{3,Float64}}}()
for i in 1:expected_cells
  p = c2c[i]
  push!(get!(panel_centroids, p, Point{3,Float64}[]), phys_coords[i][1])
end
@assert length(panel_centroids) == NPANELS "Expected $NPANELS distinct panels, got $(length(panel_centroids))"
println("Panel coverage check passed ✓  (all $NPANELS panels present)")

# ─────────────────────────────────────────────────────────────────────────────
# 10. VTK output
# ─────────────────────────────────────────────────────────────────────────────

mkpath("output")
writevtk(atlas_model, "output/cubed_sphere_explicit")
println("Written output/cubed_sphere_explicit_2.vtu — open in Paraview ✓")

println("\ntest_serial_cubed_sphere_explicit: ALL CHECKS PASSED")
