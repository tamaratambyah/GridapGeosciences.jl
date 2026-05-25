# CoarseMeshes.jl
#
# Library of canonical coarse meshes for AtlasGrid / AtlasDiscreteModel.
# Each shape is represented by a concrete subtype of AbstractCoarseMesh that
# carries the geometric parameters (radius, height, …).  Calling
# get_coarse_mesh(shape) returns a CoarseMeshInfo bundling the coarse Grid,
# per-cell local-frame corner coordinates, and default physical maps.

using GridapGeosciences
import GridapGeosciences.Fields: ForwardMap2D
import GridapGeosciences.Geometry: NPANELS, CUBE_HALF_EDGE

# ============================================================
# CoarseMeshInfo
# ============================================================

"""
    CoarseMeshInfo{Dc, G, A, M}

Bundles a coarse `Grid{Dc,Dc}` with per-cell local-frame corner coordinates
and per-chart physical maps.  Returned by `get_coarse_mesh`; consumed by
the `AtlasGrid` and `AtlasDiscreteModel` convenience constructors.

- `grid`          — coarse `UnstructuredGrid{Dc,Dc}` (topology only; coordinates are junk).
- `local_coords`  — one entry per coarse cell; `local_coords[k]` is a vector of
                    `Point{Dc}` giving the corners of chart k in its local reference frame.
- `physical_maps` — one callable per chart: `Point{Dc} → Point{Da}`.
"""
struct CoarseMeshInfo{Dc, G <: Gridap.Geometry.Grid{Dc,Dc}, A <: AbstractVector, M}
  grid          :: G
  local_coords  :: A
  physical_maps :: M

  function CoarseMeshInfo(
      grid          :: Gridap.Geometry.Grid{Dc,Dc},
      local_coords  :: A,
      physical_maps :: M,
  ) where {Dc, A <: AbstractVector, M}
    G = typeof(grid)
    new{Dc,G,A,M}(grid, local_coords, physical_maps)
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

  ref_corners  = [Point(-1.0,-1.0), Point(1.0,-1.0), Point(-1.0,1.0), Point(1.0,1.0)]
  local_coords = [ref_corners, ref_corners]

  map_C1(pt::Point{2,Float64}) =
    Point(r*cos(π*(pt[1]+1)/2), r*sin(π*(pt[1]+1)/2), h*(pt[2]+1)/2)
  map_C2(pt::Point{2,Float64}) =
    Point(r*cos(π*(pt[1]+3)/2), r*sin(π*(pt[1]+3)/2), h*(pt[2]+1)/2)

  CoarseMeshInfo(grid, local_coords, [map_C1, map_C2])
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

  panel_corners = [
    Point(-CUBE_HALF_EDGE, -CUBE_HALF_EDGE),   # BL
    Point( CUBE_HALF_EDGE, -CUBE_HALF_EDGE),   # BR
    Point(-CUBE_HALF_EDGE,  CUBE_HALF_EDGE),   # TL
    Point( CUBE_HALF_EDGE,  CUBE_HALF_EDGE),   # TR
  ]
  local_coords  = fill(panel_corners, NPANELS)
  physical_maps = [ForwardMap2D(p, m.radius) for p in 1:NPANELS]

  CoarseMeshInfo(grid, local_coords, physical_maps)
end
