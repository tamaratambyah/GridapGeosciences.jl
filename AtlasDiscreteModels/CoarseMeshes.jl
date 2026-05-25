# CoarseMeshes.jl
#
# Library of canonical coarse meshes for AtlasGrid / AtlasDiscreteModel.
# Each shape is represented by a concrete subtype of AbstractCoarseMesh that
# carries the geometric parameters (radius, height, …).  Calling
# get_coarse_mesh(shape) returns a CoarseMeshInfo bundling the coarse
# DiscreteModel (with face labels), per-cell local-frame corner coordinates,
# and default physical maps.

using GridapGeosciences
import GridapGeosciences.Fields: ForwardMap2D
import GridapGeosciences.Geometry: NPANELS, CUBE_HALF_EDGE

# ============================================================
# CoarseMeshInfo
# ============================================================

"""
    CoarseMeshInfo{Dc, Dm, A, M}

Bundles a coarse `DiscreteModel{Dc,Dc}` (with face labels) with per-cell
local-frame corner coordinates and per-chart physical maps.  Returned by
`get_coarse_mesh`; consumed by the `AtlasGrid` and `AtlasDiscreteModel`
convenience constructors.

- `model`         — coarse `DiscreteModel{Dc,Dc}` carrying topology and
                    `FaceLabeling` (node coordinates are junk — only connectivity
                    matters).  For meshes with physical boundaries (e.g. cylinder),
                    boundary edges/nodes are tagged by `get_coarse_mesh`.
- `local_coords`  — one entry per coarse cell; `local_coords[k]` is a vector of
                    `Point{Dc}` giving the corners of chart k in its local frame.
- `physical_maps` — one callable per chart: `Point{Dc} → Point{Da}`.
"""
struct CoarseMeshInfo{Dc,
                      Dm <: Gridap.Geometry.DiscreteModel{Dc,Dc},
                      A  <: AbstractVector,
                      M}
  model         :: Dm
  local_coords  :: A
  physical_maps :: M

  function CoarseMeshInfo(
      model         :: Gridap.Geometry.DiscreteModel{Dc,Dc},
      local_coords  :: A,
      physical_maps :: M,
  ) where {Dc, A <: AbstractVector, M}
    Dm = typeof(model)
    new{Dc,Dm,A,M}(model, local_coords, physical_maps)
  end
end

# ============================================================
# AbstractCoarseMesh and shape types
# ============================================================

"""
    AbstractCoarseMesh

Supertype for all canonical coarse-mesh descriptors.
Subtypes carry the geometric parameters (radius, height, …) and are passed
to `get_coarse_mesh` to obtain a `CoarseMeshInfo`.
"""
abstract type AbstractCoarseMesh end

"""
    CylinderMesh(radius=1.0, height=1.0)

Two-chart atlas for a cylinder of the given radius and height.
The two charts cover θ ∈ [0,π] (C1) and θ ∈ [π,2π] (C2), with the seam at
θ = 0 ≡ 2π encoded by periodic node sharing in the coarse mesh.
"""
struct CylinderMesh <: AbstractCoarseMesh
  radius :: Float64
  height :: Float64
  CylinderMesh(radius=1.0, height=1.0) = new(radius, height)
end

"""
    CubedSphereMesh(radius=1.0)

Six-panel atlas for a cubed sphere of the given radius.
Panels are numbered 1–6 and use the gnomonic projection `ForwardMap2D(p, radius)`,
with local (α,β) coordinates in [−π/4, π/4]².
"""
struct CubedSphereMesh <: AbstractCoarseMesh
  radius :: Float64
  CubedSphereMesh(radius=1.0) = new(radius)
end

# ============================================================
# get_coarse_mesh — CylinderMesh
# ============================================================

"""
    get_coarse_mesh(m::CylinderMesh) → CoarseMeshInfo{2}

Coarse QUAD mesh for a cylinder with radius `m.radius` and height `m.height`.

Topology (4 nodes, 2 cells). Node coordinates are junk — only connectivity
matters. Gridap Z-order per cell: BL, BR, TL, TR.

  3──4──(3)    z = height
  │  │   │
  1──2──(1)    z = 0
  [C1] [C2]

Nodes:
  1 = (0,0)   2 = (1,0)   3 = (0,1)   4 = (1,1)

Cells:
  C1 = [1, 2, 3, 4]   BL=1 BR=2 TL=3 TR=4   θ ∈ [0,  π]
  C2 = [2, 1, 4, 3]   BL=2 BR=1 TL=4 TR=3   θ ∈ [π, 2π]

Node sharing: C2's BR=node1=C1's BL and C2's TR=node3=C1's TL.
This encodes the periodic identification θ=0 ≡ θ=2π (the seam).

Local frame for both charts: (s,t) ∈ [−1, 1]².
Physical maps (s,t) → (X,Y,Z):
  C1: (r·cos(π(s+1)/2), r·sin(π(s+1)/2), h(t+1)/2)   r = radius, h = height
  C2: (r·cos(π(s+3)/2), r·sin(π(s+3)/2), h(t+1)/2)

Face labels (entity 1 = "interior", entity 2 = "boundary" reserved by Gridap):
  "bottom" = entity 3: edge (1,2) and nodes 1, 2  (z = 0 circle)
  "top"    = entity 4: edge (3,4) and nodes 3, 4  (z = height circle)
"""
function get_coarse_mesh(m::CylinderMesh)
  r, h = m.radius, m.height

  node_coords = Vector{Point{2,Float64}}([
    Point(0.0, 0.0),   # 1
    Point(1.0, 0.0),   # 2
    Point(0.0, 1.0),   # 3
    Point(1.0, 1.0),   # 4
  ])

  cell_node_data = Int32[1,2,3,4,  2,1,4,3]
  cell_node_ptrs = Int32[1,5,9]
  cell_node_ids  = Gridap.Arrays.Table(cell_node_data, cell_node_ptrs)
  cell_types     = Int32[1,1]
  reffe          = Gridap.ReferenceFEs.LagrangianRefFE(Float64, QUAD, 1)
  grid = Gridap.Geometry.UnstructuredGrid(
    node_coords, cell_node_ids, [reffe], cell_types, Gridap.Geometry.NonOriented())

  # ── Topology + face labels with "bottom" and "top" tags ───────────────────
  # All 4 edges are topologically interior (periodic mesh); get_isboundary_face
  # returns false for all.  The z=0 and z=height circles are physical boundaries
  # tagged explicitly.  FaceLabeling(topo) uses entity 1 for all interior faces;
  # we override with entity 3 (bottom) and 4 (top).  Entities 1 and 2 are
  # reserved by Gridap ("interior" and "boundary").
  #
  # Edges are identified by their vertex-pair rather than by positional index
  # so the labeling is robust to Gridap changing its edge-enumeration order.
  # Bottom circle: edge connecting nodes 1 and 2 (z = 0).
  # Top circle:    edge connecting nodes 3 and 4 (z = height).
  topo   = Gridap.Geometry.UnstructuredGridTopology(grid)
  labels = Gridap.Geometry.FaceLabeling(topo)   # all interior entities = 1

  edge_to_vert = Gridap.Geometry.get_faces(topo, 1, 0)
  n_edges      = Gridap.Geometry.num_faces(topo, 1)
  for e in 1:n_edges
    vs = sort(collect(Int32, edge_to_vert[e]))
    if vs == Int32[1, 2]
      labels.d_to_dface_to_entity[2][e] = Int32(3)   # bottom circle
    elseif vs == Int32[3, 4]
      labels.d_to_dface_to_entity[2][e] = Int32(4)   # top circle
    end
  end

  # Propagate edge entities to their endpoint nodes.
  for e in 1:n_edges
    entity = labels.d_to_dface_to_entity[2][e]
    if entity == Int32(3) || entity == Int32(4)
      for v in edge_to_vert[e]
        labels.d_to_dface_to_entity[1][v] = entity
      end
    end
  end

  Gridap.Geometry.add_tag!(labels, "bottom", [3])
  Gridap.Geometry.add_tag!(labels, "top",    [4])

  model = Gridap.Geometry.UnstructuredDiscreteModel(grid, topo, labels)

  ref_corners  = [Point(-1.0,-1.0), Point(1.0,-1.0), Point(-1.0,1.0), Point(1.0,1.0)]
  local_coords = [ref_corners, ref_corners]

  map_C1(pt::Point{2,Float64}) =
    Point(r*cos(π*(pt[1]+1)/2), r*sin(π*(pt[1]+1)/2), h*(pt[2]+1)/2)
  map_C2(pt::Point{2,Float64}) =
    Point(r*cos(π*(pt[1]+3)/2), r*sin(π*(pt[1]+3)/2), h*(pt[2]+1)/2)

  CoarseMeshInfo(model, local_coords, [map_C1, map_C2])
end

# ============================================================
# get_coarse_mesh — CubedSphereMesh
# ============================================================

"""
    get_coarse_mesh(m::CubedSphereMesh) → CoarseMeshInfo{2}

Coarse QUAD mesh for a cubed sphere with radius `m.radius`.

Topology (8 nodes, 6 cells). Node coordinates are junk — only connectivity
matters. Gridap Z-order per cell: BL, BR, TL, TR.

       x=−2                              x=2
  y=2   5──────────────────────────────────6   C2
        │╲              C4               ╱│
  y=1   │  8──────────────────────────7  │
        │  │                          │  │
        │C6│         C5               │C3│
        │  │                          │  │
  y=−1  │  1──────────────────────────2  │
        │╱              C1               ╲│
  y=−2  3──────────────────────────────────4

Nodes:
  1=(−1,−1)  2=(1,−1)  3=(−2,−2)  4=(2,−2)
  5=(−2, 2)  6=(2, 2)  7=(1, 1)   8=(−1, 1)

Cells (same connectivity as _CCAM_panel_wise_node_ids):
  C1 = [1,2,3,4]   BL=1 BR=2 TL=3 TR=4
  C2 = [3,4,5,6]   BL=3 BR=4 TL=5 TR=6
  C3 = [2,7,4,6]   BL=2 BR=7 TL=4 TR=6
  C4 = [8,5,7,6]   BL=8 BR=5 TL=7 TR=6
  C5 = [1,8,2,7]   BL=1 BR=8 TL=2 TR=7
  C6 = [1,3,8,5]   BL=1 BR=3 TL=8 TR=5

The 8 nodes correspond to the 8 corners of the cube. Node sharing encodes the
12 shared edges between the 6 faces.

Local frame for all charts: (α,β) ∈ [−π/4, π/4]².
Physical maps: `ForwardMap2D(p, radius)` for p = 1 … 6 (gnomonic projection).

Face labels: the cubed sphere is a closed manifold — no topological boundary.
FaceLabeling(topo) assigns entity 1 ("interior") to all faces, which also
satisfies the p4est precondition that all face entity ids be positive.
# NOTE: The legacy serial function _create_parametric_octree_dmodel_coarse_model
# explicitly calls add_tag!(labels, "boundary", [1]).  That call is UNNECESSARY
# (and would fail with current Gridap because "boundary" already exists in the
# default FaceLabeling).  The p4est precondition is entity id > 0, which is
# already satisfied here.
"""
function get_coarse_mesh(m::CubedSphereMesh)
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

  cell_node_data = Int32[1,2,3,4, 3,4,5,6, 2,7,4,6, 8,5,7,6, 1,8,2,7, 1,3,8,5]
  cell_node_ptrs = Int32[1,5,9,13,17,21,25]
  cell_node_ids  = Gridap.Arrays.Table(cell_node_data, cell_node_ptrs)
  cell_types     = Int32[1,1,1,1,1,1]
  reffe          = Gridap.ReferenceFEs.LagrangianRefFE(Float64, QUAD, 1)
  grid = Gridap.Geometry.UnstructuredGrid(
    node_coords, cell_node_ids, [reffe], cell_types, Gridap.Geometry.NonOriented())

  topo   = Gridap.Geometry.UnstructuredGridTopology(grid)
  labels = Gridap.Geometry.FaceLabeling(topo)   # all entities = 1 (closed manifold)
  model  = Gridap.Geometry.UnstructuredDiscreteModel(grid, topo, labels)

  panel_corners = [
    Point(-CUBE_HALF_EDGE, -CUBE_HALF_EDGE),   # BL
    Point( CUBE_HALF_EDGE, -CUBE_HALF_EDGE),   # BR
    Point(-CUBE_HALF_EDGE,  CUBE_HALF_EDGE),   # TL
    Point( CUBE_HALF_EDGE,  CUBE_HALF_EDGE),   # TR
  ]
  local_coords  = fill(panel_corners, NPANELS)
  physical_maps = [ForwardMap2D(p, m.radius) for p in 1:NPANELS]

  CoarseMeshInfo(model, local_coords, physical_maps)
end
