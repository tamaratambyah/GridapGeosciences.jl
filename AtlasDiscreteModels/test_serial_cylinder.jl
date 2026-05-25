# test_serial_cylinder.jl
#
# Tests AtlasGrid / AtlasDiscreteModel on a 2-chart cylinder atlas.
#
# Coarse mesh: 2 QUAD cells with 4 nodes (not 6). The right boundary of C2
# is topologically identified with the left boundary of C1, forming the seam
# of the cylinder.  Gridap Z-order per cell: BL, BR, TL, TR.
#
#   3──4──(3)    (z = 1)
#   │  │   │
#   1──2──(1)    (z = 0)
#   [C1] [C2]
#
# Parenthesised nodes are the same objects as nodes 1 and 3 — node sharing
# encodes the periodic identification θ = 0 ≡ θ = 2π.
#
# C1: nodes (1, 2, 3, 4) — BL=1, BR=2, TL=3, TR=4   θ ∈ [0, π]
# C2: nodes (2, 1, 4, 3) — BL=2, BR=1, TL=4, TR=3   θ ∈ [π, 2π]
#      C2's BR=node1 = C1's BL and C2's TR=node3 = C1's TL: the seam.
#
# Each chart has its own local coordinate system: standard QUAD ref element [-1,1]².
#
# Physical maps (local [-1,1]² → cylinder surface ℝ³):
#   map_C1 : (s,t) → (cos(π(s+1)/2), sin(π(s+1)/2), (t+1)/2)   θ ∈ [0,π]
#   map_C2 : (s,t) → (cos(π(s+3)/2), sin(π(s+3)/2), (t+1)/2)   θ ∈ [π,2π]
#
# Run:
#   julia --project=. AtlasDiscreteModels/test_serial_cylinder.jl

using Gridap
using Gridap.Geometry, Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Adaptivity, Gridap.Helpers
using FillArrays

include("AtlasDiscreteModels.jl")

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Build coarse 2-cell QUAD mesh with periodic identification at the seam
# ─────────────────────────────────────────────────────────────────────────────

node_coords = Vector{Point{2,Float64}}([
  Point(0.0, 0.0),   # 1
  Point(1.0, 0.0),   # 2
  Point(0.0, 1.0),   # 3
  Point(1.0, 1.0),   # 4
])

# C1: (BL=1, BR=2, TL=3, TR=4)
# C2: (BL=2, BR=1, TL=4, TR=3) — BR and TR wrap back to C1's left edge
cell_node_data = Int32[1,2,3,4,  2,1,4,3]
cell_node_ptrs = Int32[1,5,9]
cell_node_ids  = Table(cell_node_data, cell_node_ptrs)
cell_types     = Int32[1,1]
reffe          = LagrangianRefFE(Float64, QUAD, 1)
coarse_grid    = UnstructuredGrid(node_coords, cell_node_ids, [reffe], cell_types, NonOriented())
coarse_topo    = UnstructuredGridTopology(node_coords, cell_node_ids, cell_types, [QUAD], NonOriented())
coarse_labels  = FaceLabeling(coarse_topo)
coarse_model   = UnstructuredDiscreteModel(coarse_grid, coarse_topo, coarse_labels)

println("Coarse mesh: $(num_cells(coarse_model)) cells, $(num_nodes(coarse_model)) nodes")
@assert num_nodes(coarse_model) == 4 "Expected 4 nodes (periodic mesh), got $(num_nodes(coarse_model))"

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Local coordinates for each coarse chart
#     Both charts use the standard QUAD reference element corners [-1,1]².
#     Gridap Z-order: BL, BR, TL, TR
# ─────────────────────────────────────────────────────────────────────────────

ref_corners = [Point(-1.0,-1.0), Point(1.0,-1.0), Point(-1.0,1.0), Point(1.0,1.0)]
coarse_local_coords = [ref_corners, ref_corners]

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Physical maps: local [-1,1]² → cylinder surface
# ─────────────────────────────────────────────────────────────────────────────

map_C1(pt::Point{2,Float64}) = Point(cos(π*(pt[1]+1)/2), sin(π*(pt[1]+1)/2), (pt[2]+1)/2)
map_C2(pt::Point{2,Float64}) = Point(cos(π*(pt[1]+3)/2), sin(π*(pt[1]+3)/2), (pt[2]+1)/2)

physical_maps = [map_C1, map_C2]

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Build AtlasDiscreteModel
# ─────────────────────────────────────────────────────────────────────────────

atlas_model = AtlasDiscreteModel(coarse_model, coarse_local_coords, physical_maps, 3)
atlas_grid  = get_atlas_grid(atlas_model)

println("\nAtlasGrid:  Dc=$(num_cell_dims(atlas_grid)), Da=$(get_ambient_dim(atlas_grid))")
println("            cells = $(Gridap.Geometry.num_cells(atlas_grid))")
println("AtlasModel: built ✓")

# ─────────────────────────────────────────────────────────────────────────────
# 5.  Sanity checks on structure
# ─────────────────────────────────────────────────────────────────────────────

@assert Gridap.Geometry.num_cells(atlas_grid) == 128   # 2 coarse × 4^3

c2c = get_cell_to_chart(atlas_grid)
println("cell_to_chart (first 8): ", c2c[1:8])
@assert all(c2c[1:64]   .== 1)
@assert all(c2c[65:128] .== 2)

# get_cell_coordinates returns LOCAL 2D coords
local_coords = Gridap.Geometry.get_cell_coordinates(atlas_grid)
@assert length(local_coords) == 128

println("Local coords cell 1:  ", local_coords[1])
println("Local coords cell 65: ", local_coords[65])

# All local coords must lie within [-1,1]²
for i in 1:Gridap.Geometry.num_cells(atlas_grid)
  for pt in local_coords[i]
    @assert -1.0 - 1e-12 ≤ pt[1] ≤ 1.0 + 1e-12 "cell $i: local x=$(pt[1]) out of [-1,1]"
    @assert -1.0 - 1e-12 ≤ pt[2] ≤ 1.0 + 1e-12 "cell $i: local y=$(pt[2]) out of [-1,1]"
  end
end
println("Local coord range check passed ✓")

# ─────────────────────────────────────────────────────────────────────────────
# 6.  Sanity check on physical coords (via _local_to_physical)
# ─────────────────────────────────────────────────────────────────────────────

phys_coords = _local_to_physical(
  atlas_grid.cell_local_coords,
  atlas_grid.cell_to_chart,
  atlas_grid.physical_maps,
)
@assert length(phys_coords) == 128

for i in 1:Gridap.Geometry.num_cells(atlas_grid)
  for pt in phys_coords[i]
    r = sqrt(pt[1]^2 + pt[2]^2)
    @assert isapprox(r, 1.0; atol=1e-12) "cell $i: r=$r ≠ 1  pt=$pt"
    @assert 0.0 - 1e-12 ≤ pt[3] ≤ 1.0 + 1e-12 "cell $i: z=$(pt[3]) out of [0,1]"
  end
end
println("Physical coord cylinder check passed ✓")

# ─────────────────────────────────────────────────────────────────────────────
# 7.  Seam check: map_C1(s=1, t) == map_C2(s=-1, t) for all t (θ=π is shared)
#     and map_C1(s=-1, t) == map_C2(s=1, t) for all t (θ=0 ≡ θ=2π, the seam)
# ─────────────────────────────────────────────────────────────────────────────

for t in range(-1.0, 1.0; length=9)
  interface = map_C1(Point(1.0, t))   # θ = π, right edge of C1
  @assert interface ≈ map_C2(Point(-1.0, t)) "Interface mismatch at t=$t"
  seam_C1 = map_C1(Point(-1.0, t))   # θ = 0
  seam_C2 = map_C2(Point(1.0, t))    # θ = 2π ≡ 0
  @assert seam_C1 ≈ seam_C2 "Seam mismatch at t=$t: C1=$(seam_C1) C2=$(seam_C2)"
end
println("Seam continuity check passed ✓")

# ─────────────────────────────────────────────────────────────────────────────
# 8.  Write VTK
# ─────────────────────────────────────────────────────────────────────────────

mkpath("output")
writevtk(atlas_model, "output/cylinder")
println("Written output/cylinder_2.vtu — open in Paraview ✓")
