# test_serial_cylinder.jl
#
# Tests AtlasGrid / AtlasDiscreteModel on a 2-chart cylinder atlas using the
# CylinderMesh coarse mesh library.
#
# Coarse mesh (from CoarseMeshes.jl): 4 nodes, 2 QUAD cells with periodic
# identification at the seam.  See get_coarse_mesh(::CylinderMesh) for the
# full topology diagram.
#
# Run:
#   julia --project=. AtlasDiscreteModels/test_serial_cylinder.jl

using Gridap
using GridapGeosciences

include("AtlasDiscreteModels.jl")

const RADIUS = 1.0
const HEIGHT = 1.0
const NUM_REF = 3   # 2 coarse × 4^3 = 128 fine cells

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Build AtlasDiscreteModel via CylinderMesh
# ─────────────────────────────────────────────────────────────────────────────

atlas_model = AtlasDiscreteModel(CylinderMesh(RADIUS, HEIGHT), NUM_REF)
atlas_grid  = get_atlas_grid(atlas_model)

println("AtlasGrid:  Dc=$(Gridap.Geometry.num_cell_dims(atlas_grid)), Da=$(get_ambient_dim(atlas_grid))")
println("            cells = $(Gridap.Geometry.num_cells(atlas_grid))")
println("AtlasModel: built ✓")

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Structural checks
# ─────────────────────────────────────────────────────────────────────────────

@assert Gridap.Geometry.num_cells(atlas_grid) == 2 * 4^NUM_REF

# get_cell_coordinates returns LOCAL 2D coords in [-1,1]²
local_coords = Gridap.Geometry.get_cell_coordinates(atlas_grid)
@assert length(local_coords) == 2 * 4^NUM_REF
for i in 1:Gridap.Geometry.num_cells(atlas_grid)
  for pt in local_coords[i]
    @assert -1.0 - 1e-12 ≤ pt[1] ≤ 1.0 + 1e-12 "cell $i: local x=$(pt[1]) out of [-1,1]"
    @assert -1.0 - 1e-12 ≤ pt[2] ≤ 1.0 + 1e-12 "cell $i: local y=$(pt[2]) out of [-1,1]"
  end
end
println("Local coord range check passed ✓")

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Physical coords: r = RADIUS, z ∈ [0, HEIGHT]
# ─────────────────────────────────────────────────────────────────────────────

phys_coords = _local_to_physical(
  atlas_grid.cell_local_coords,
  atlas_grid.cell_physical_maps,
)
for i in 1:Gridap.Geometry.num_cells(atlas_grid)
  for pt in phys_coords[i]
    r = sqrt(pt[1]^2 + pt[2]^2)
    @assert isapprox(r, RADIUS; atol=1e-12) "cell $i: r=$r ≠ $RADIUS"
    @assert 0.0 - 1e-12 ≤ pt[3] ≤ HEIGHT + 1e-12 "cell $i: z=$(pt[3]) out of [0,$HEIGHT]"
  end
end
println("Physical coord cylinder check passed ✓")

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Seam continuity via the maps stored in CoarseMeshInfo
# ─────────────────────────────────────────────────────────────────────────────

map_C1, map_C2 = get_coarse_mesh(CylinderMesh(RADIUS, HEIGHT)).physical_maps
for t in range(-1.0, 1.0; length=9)
  @assert map_C1(Point(1.0, t))  ≈ map_C2(Point(-1.0, t)) "Interface mismatch at t=$t"
  @assert map_C1(Point(-1.0, t)) ≈ map_C2(Point(1.0, t))  "Seam mismatch at t=$t"
end
println("Seam continuity check passed ✓")

# ─────────────────────────────────────────────────────────────────────────────
# 5.  Face label propagation: "bottom" and "top" survive refinement
# ─────────────────────────────────────────────────────────────────────────────

fine_labels   = Gridap.Geometry.get_face_labeling(atlas_model)
edge_entities = Gridap.Geometry.get_face_entity(fine_labels, 1)
bottom_tag    = Gridap.Geometry.get_tag_from_name(fine_labels, "bottom")
top_tag       = Gridap.Geometry.get_tag_from_name(fine_labels, "top")
bottom_entity = fine_labels.tag_to_entities[bottom_tag][1]
top_entity    = fine_labels.tag_to_entities[top_tag][1]
bottom_edges  = findall(==(bottom_entity), edge_entities)
top_edges     = findall(==(top_entity),    edge_entities)
@assert !isempty(bottom_edges) "\"bottom\" tag not found on fine mesh edges"
@assert !isempty(top_edges)    "\"top\" tag not found on fine mesh edges"
# Each coarse bottom/top edge refines into 2^NUM_REF fine edges.
@assert length(bottom_edges) == 2^NUM_REF "Expected $(2^NUM_REF) bottom edges, got $(length(bottom_edges))"
@assert length(top_edges)    == 2^NUM_REF "Expected $(2^NUM_REF) top edges, got $(length(top_edges))"
println("Face label propagation check passed ✓  ($(2^NUM_REF) bottom, $(2^NUM_REF) top edges)")

# ─────────────────────────────────────────────────────────────────────────────
# 6.  VTK output
# ─────────────────────────────────────────────────────────────────────────────

mkpath("output")
writevtk(atlas_model, "output/cylinder")
println("Written output/cylinder_2.vtu — open in Paraview ✓")
