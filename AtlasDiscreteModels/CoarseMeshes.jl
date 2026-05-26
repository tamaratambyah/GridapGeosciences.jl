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
    CoarseMeshInfo{Dc, Dm, A, M, G}

Bundles a coarse `DiscreteModel{Dc,Dc}` (with face labels) with per-cell
local-frame corner coordinates, per-chart physical maps, and per-chart metric fields.
Returned by `get_coarse_mesh`; consumed by the `AtlasGrid` and `AtlasDiscreteModel`
convenience constructors.

- `model`         — coarse `DiscreteModel{Dc,Dc}` carrying topology and
                    `FaceLabeling` (node coordinates are junk — only connectivity
                    matters).  For meshes with physical boundaries (e.g. cylinder),
                    boundary edges/nodes are tagged by `get_coarse_mesh`.
- `local_coords`  — one entry per coarse cell; `local_coords[k]` is a vector of
                    `Point{Dc}` giving the corners of chart k in its local frame.
- `physical_maps` — one `Field` per chart: `Point{Dc} → Point{Da}`.
- `metric_fields` — one `Field` per chart: `Point{Dc} → SymTensorValue{Dc}`,
                    the pullback metric `g = JᵀJ`.  Computed lazily from
                    `physical_maps` via `Operation(JtJ())(gradient(φ))` when not
                    provided explicitly as an analytic expression.
"""
struct CoarseMeshInfo{Dc,
                      Dm <: Gridap.Geometry.DiscreteModel{Dc,Dc},
                      A  <: AbstractVector,
                      M,
                      G}
  model         :: Dm
  local_coords  :: A
  physical_maps :: M
  metric_fields :: G

  function CoarseMeshInfo(
      model         :: Gridap.Geometry.DiscreteModel{Dc,Dc},
      local_coords  :: A,
      physical_maps :: M,
      metric_fields :: G,
  ) where {Dc, A <: AbstractVector, M, G}
    Dm = typeof(model)
    new{Dc,Dm,A,M,G}(model, local_coords, physical_maps, metric_fields)
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

"""
    MobiusStripMesh(radius=1.0, half_width=0.3)

Two-chart atlas for a Möbius strip with major radius `radius` and half-width `half_width`.
The two charts cover θ ∈ [0,π] (C1) and θ ∈ [π,2π] (C2). The half-twist is encoded
topologically: the right edge of C2 is identified with the left edge of C1 reversed.
"""
struct MobiusStripMesh <: AbstractCoarseMesh
  radius     :: Float64
  half_width :: Float64
  MobiusStripMesh(radius=1.0, half_width=0.3) = new(radius, half_width)
end

# ============================================================
# Concrete Field subtypes for physical embeddings
# ============================================================
#
# Plain Julia closures create a Vector{Function} (abstract element type) when
# collected into an array, causing dynamic dispatch through all downstream
# lazy_map chains.  Concrete Field subtypes give a fully typed array, keeping
# Reindex / lazy_map type-stable.
#
# Each type must implement the full Gridap Field interface:
#   evaluate!(cache, f, x::VectorValue{2})           — single point
#   evaluate!(cache, f, xs::AbstractArray{...})       — vector of points
#   return_cache(f, xs::AbstractArray{...})           — preallocate output
# plus the FieldGradient methods for the Jacobian (needed by integration).
#
# Gridap gradient convention (see github.com/gridap/Gridap.jl/issues/822):
# FieldGradient{1,F} evaluated at x returns ∇F(x) = JT, stored as
# TensorValue{Dc,Da,Float64} with column-major data (dX/ds,dX/dt, dY/ds,dY/dt, dZ/ds,dZ/dt).

"""
    CylinderChartMap(radius, height, theta_offset)

Gridap `Field` for one chart of the cylinder atlas.
`theta_offset = 1.0` → θ ∈ [0,π] (C1);  `3.0` → θ ∈ [π,2π] (C2).
Maps (s,t) ∈ [−1,1]² to (r·cos θ, r·sin θ, h(t+1)/2) with θ = π(s+offset)/2.
"""
struct CylinderChartMap <: Field
  radius       :: Float64
  height       :: Float64
  theta_offset :: Float64
end

function Gridap.Arrays.evaluate!(cache, m::CylinderChartMap, x::VectorValue{2})
  θ = π*(x[1] + m.theta_offset)/2
  Point(m.radius*cos(θ), m.radius*sin(θ), m.height*(x[2]+1)/2)
end

function Gridap.Arrays.return_cache(m::CylinderChartMap, xs::AbstractArray{<:VectorValue{2}})
  similar(xs, VectorValue{3,Float64})
end
function Gridap.Arrays.evaluate!(cache, m::CylinderChartMap, xs::AbstractArray{<:VectorValue{2}})
  cache .= evaluate!.(nothing, Ref(m), xs)
  cache
end

# Jacobian transpose JT = TensorValue{2,3}: column-major (dX/ds,dX/dt, dY/ds,dY/dt, dZ/ds,dZ/dt)
function _cylinder_jacobian_T(m::CylinderChartMap, x::VectorValue{2})
  θ = π*(x[1] + m.theta_offset)/2
  TensorValue{2,3,Float64}(
    -m.radius*sin(θ)*(π/2), 0.0,
     m.radius*cos(θ)*(π/2), 0.0,
     0.0,                   m.height/2,
  )
end

function Gridap.Arrays.return_cache(
    f::FieldGradient{1,<:CylinderChartMap}, xs::AbstractArray{<:VectorValue{2}})
  CachedArray(similar(xs, TensorValue{2,3,Float64}))
end
function Gridap.Arrays.evaluate!(
    cache, f::FieldGradient{1,<:CylinderChartMap}, xs::AbstractArray{<:VectorValue{2}})
  setsize!(cache, size(xs))
  cache.array .= _cylinder_jacobian_T.(Ref(f.object), xs)
  cache.array
end
function Gridap.Arrays.evaluate!(
    cache, f::FieldGradient{1,<:CylinderChartMap}, x::VectorValue{2})
  _cylinder_jacobian_T(f.object, x)
end

"""
    MobiusChartMap(radius, half_width, theta_offset)

Gridap `Field` for one chart of the Möbius strip atlas.
`theta_offset = 1.0` → θ ∈ [0,π] (C1);  `3.0` → θ ∈ [π,2π] (C2).
Maps (s,t) ∈ [−1,1]² to the standard Möbius embedding with half-angle twist.
"""
struct MobiusChartMap <: Field
  radius       :: Float64
  half_width   :: Float64
  theta_offset :: Float64
end

function Gridap.Arrays.evaluate!(cache, m::MobiusChartMap, x::VectorValue{2})
  θ = π*(x[1] + m.theta_offset)/2
  ρ = m.radius + m.half_width*x[2]*cos(θ/2)
  Point(ρ*cos(θ), ρ*sin(θ), m.half_width*x[2]*sin(θ/2))
end

function Gridap.Arrays.return_cache(m::MobiusChartMap, xs::AbstractArray{<:VectorValue{2}})
  similar(xs, VectorValue{3,Float64})
end
function Gridap.Arrays.evaluate!(cache, m::MobiusChartMap, xs::AbstractArray{<:VectorValue{2}})
  cache .= evaluate!.(nothing, Ref(m), xs)
  cache
end

# Jacobian transpose JT = TensorValue{2,3}: column-major (dX/ds,dX/dt, dY/ds,dY/dt, dZ/ds,dZ/dt)
function _mobius_jacobian_T(m::MobiusChartMap, x::VectorValue{2})
  θ = π*(x[1] + m.theta_offset)/2
  t = x[2]
  W = m.half_width
  ρ = m.radius + W*t*cos(θ/2)
  dρds = -W*t*sin(θ/2)*(π/4)
  dρdt =  W*cos(θ/2)
  TensorValue{2,3,Float64}(
    dρds*cos(θ) - ρ*sin(θ)*(π/2),  dρdt*cos(θ),    # dX/ds, dX/dt
    dρds*sin(θ) + ρ*cos(θ)*(π/2),  dρdt*sin(θ),    # dY/ds, dY/dt
    W*t*cos(θ/2)*(π/4),             W*sin(θ/2),     # dZ/ds, dZ/dt
  )
end

function Gridap.Arrays.return_cache(
    f::FieldGradient{1,<:MobiusChartMap}, xs::AbstractArray{<:VectorValue{2}})
  CachedArray(similar(xs, TensorValue{2,3,Float64}))
end
function Gridap.Arrays.evaluate!(
    cache, f::FieldGradient{1,<:MobiusChartMap}, xs::AbstractArray{<:VectorValue{2}})
  setsize!(cache, size(xs))
  cache.array .= _mobius_jacobian_T.(Ref(f.object), xs)
  cache.array
end
function Gridap.Arrays.evaluate!(
    cache, f::FieldGradient{1,<:MobiusChartMap}, x::VectorValue{2})
  _mobius_jacobian_T(f.object, x)
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

  # Coordinate values are unused — AtlasGrid replaces them with local_coords.
  # Any distinct values that give a valid non-degenerate mesh work here.
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

  ref_corners   = [Point(-1.0,-1.0), Point(1.0,-1.0), Point(-1.0,1.0), Point(1.0,1.0)]
  local_coords  = [ref_corners, ref_corners]
  physical_maps = [CylinderChartMap(r, h, 1.0), CylinderChartMap(r, h, 3.0)]
  metric_fields = [Operation(JtJ())(gradient(φ)) for φ in physical_maps]

  CoarseMeshInfo(model, local_coords, physical_maps, metric_fields)
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
  # Coordinate values are unused — AtlasGrid replaces them with local_coords.
  # Any distinct values that give a valid non-degenerate mesh work here.
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
  metric_fields = [Operation(JtJ())(gradient(φ)) for φ in physical_maps]

  CoarseMeshInfo(model, local_coords, physical_maps, metric_fields)
end

# ============================================================
# get_coarse_mesh — MobiusStripMesh
# ============================================================

"""
    get_coarse_mesh(m::MobiusStripMesh) → CoarseMeshInfo{2}

Coarse QUAD mesh for a Möbius strip with major radius `m.radius` and half-width
`m.half_width`.

Topology (4 nodes, 2 cells). Node coordinates are junk — only connectivity matters.
Gridap Z-order per cell: BL, BR, TL, TR.

  3 - 4 - 1
  |   |   |
  1 - 2 - 3
  [C1] [C2]

Cells:
  C1 = [1,2,3,4]   BL=1 BR=2 TL=3 TR=4   θ ∈ [0,  π]
  C2 = [2,3,4,1]   BL=2 BR=3 TL=4 TR=1   θ ∈ [π, 2π]

The left edge of C1 {1,3} is identified with the right edge of C2 {3,1}: same node
set, reversed orientation — this encodes the half-twist.  Compare with the cylinder
(C2=[2,1,4,3]) where the seam nodes share the same height.

Local frame for both charts: (s,t) ∈ [−1,1]², s = angular direction, t = width.
Physical maps (s,t) → (X,Y,Z), with R = major radius, W = half_width:
  C1: θ = π(s+1)/2 ∈ [0,π]
  C2: θ = π(s+3)/2 ∈ [π,2π]
  both: ((R + W·t·cos(θ/2))·cos(θ),  (R + W·t·cos(θ/2))·sin(θ),  W·t·sin(θ/2))

Seam continuity:
  Interior (θ = π):   map_C1(1, t)  = map_C2(−1, t)
  Twist    (θ = 0≡2π): map_C1(−1, t) = map_C2(1, −t)   [t ↦ −t encodes the half-twist]
"""
function get_coarse_mesh(m::MobiusStripMesh)
  R, W = m.radius, m.half_width

  # Coordinate values are unused — AtlasGrid replaces them with local_coords.
  # Any distinct values that give a valid non-degenerate mesh work here.
  node_coords = Vector{Point{2,Float64}}([
    Point(0.0, 0.0),   # 1
    Point(1.0, 0.0),   # 2
    Point(0.0, 1.0),   # 3
    Point(1.0, 1.0),   # 4
  ])

  cell_node_data = Int32[1,2,3,4,  2,3,4,1]
  cell_node_ptrs = Int32[1,5,9]
  cell_node_ids  = Gridap.Arrays.Table(cell_node_data, cell_node_ptrs)
  cell_types     = Int32[1,1]
  reffe          = Gridap.ReferenceFEs.LagrangianRefFE(Float64, QUAD, 1)
  grid = Gridap.Geometry.UnstructuredGrid(
    node_coords, cell_node_ids, [reffe], cell_types, Gridap.Geometry.NonOriented())

  topo   = Gridap.Geometry.UnstructuredGridTopology(grid)
  labels = Gridap.Geometry.FaceLabeling(topo)
  model  = Gridap.Geometry.UnstructuredDiscreteModel(grid, topo, labels)

  ref_corners   = [Point(-1.0,-1.0), Point(1.0,-1.0), Point(-1.0,1.0), Point(1.0,1.0)]
  local_coords  = [ref_corners, ref_corners]
  physical_maps = [MobiusChartMap(R, W, 1.0), MobiusChartMap(R, W, 3.0)]
  metric_fields = [Operation(JtJ())(gradient(φ)) for φ in physical_maps]

  CoarseMeshInfo(model, local_coords, physical_maps, metric_fields)
end
