# SerialAtlasModels.jl
#
# Serial (no-MPI) definitions of AtlasGrid, AtlasDiscreteModel, and the
# _transform_alpha_beta_to_3d helper.
#
# Copied verbatim from SBAtlasModels.jl (lines 22–243) with all parallel
# imports (GridapDistributed, GridapP4est, PartitionedArrays, MPI) removed.
# The parallel-only AtlasOctreeDistributedDiscreteModel and its constructor
# remain in SBAtlasModels.jl only.

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Helpers
using FillArrays

# ============================================================
# AtlasGrid
# ============================================================

"""
    AtlasGrid{Dc,Da,Tp,G,A,M,O} <: Gridap.Geometry.Grid{Dc,Da}

A DG-style grid representing a `Dc`-dimensional manifold embedded in a `Da`-dimensional
ambient space (e.g. Dc=2, Da=3 for a 2D sphere surface in ℝ³).

# Fields
- `topo_grid::G`         — Purely topological grid in `Dc`-dimensional parameter space
                           (coordinates are junk; only connectivity is used).
- `cell_coordinates::A`  — Per-cell ambient coordinates (Da-dimensional). DG-style:
                           each cell has its own vector of `Da`-dim corner points.
                           No shared nodes between cells.
- `cell_to_chart::Vector{Int}` — For each fine cell, the index of its coarse chart panel.
- `physical_maps::M`     — One chart map (Dc → Da) per coarse panel, used to produce
                           `cell_coordinates` from the parametric quadrant positions.
- `reffe`                — The single `LagrangianRefFE{Dc}` shared by all cells.
- `orientation_style::O` — `Gridap.Geometry.Oriented()` or `NonOriented()`.

# Notes
- `get_node_coordinates` is **not implemented**: there are no globally shared nodes.
  Use `get_cell_coordinates` to access the Da-dimensional corner positions per cell.
- `get_cell_map` is automatically derived from `get_cell_coordinates` (via the
  Gridap default implementation), returning Da-dimensional cell maps suitable for
  FEM integration on the manifold.
"""
struct AtlasGrid{Dc, Da, Tp,
                 G <: Gridap.Geometry.Grid,
                 A <: AbstractVector,
                 M,
                 O <: Gridap.Geometry.OrientationStyle} <: Gridap.Geometry.Grid{Dc,Da}
  topo_grid        :: G
  cell_coordinates :: A
  cell_to_chart    :: Vector{Int}
  physical_maps    :: M
  reffe            :: Gridap.ReferenceFEs.LagrangianRefFE{Dc}
  orientation_style :: O

  function AtlasGrid(
    topo_grid        :: Gridap.Geometry.Grid{Dc,Dc},
    cell_coordinates :: AbstractVector{<:AbstractVector{<:Point{Da,Tp}}},
    cell_to_chart    :: AbstractVector{<:Integer},
    physical_maps,
    reffe            :: Gridap.ReferenceFEs.LagrangianRefFE{Dc},
    orientation_style :: Gridap.Geometry.OrientationStyle
  ) where {Dc,Da,Tp}
    n = Gridap.Geometry.num_cells(topo_grid)
    @check length(cell_coordinates) == n
        "cell_coordinates length ($(length(cell_coordinates))) ≠ num_cells ($(n))"
    @check length(cell_to_chart) == n
        "cell_to_chart length ($(length(cell_to_chart))) ≠ num_cells ($(n))"
    G = typeof(topo_grid)
    A = typeof(cell_coordinates)
    M = typeof(physical_maps)
    O = typeof(orientation_style)
    new{Dc,Da,Tp,G,A,M,O}(
      topo_grid,
      cell_coordinates,
      Vector{Int}(cell_to_chart),
      physical_maps,
      reffe,
      orientation_style,
    )
  end
end

# ----------------------------------------------------------
# Gridap.Geometry.Grid{Dc,Da} interface
# ----------------------------------------------------------

Gridap.Geometry.OrientationStyle(
  ::Type{<:AtlasGrid{Dc,Da,Tp,G,A,M,O}}) where {Dc,Da,Tp,G,A,M,O} = O()

Gridap.Geometry.num_cells(g::AtlasGrid) = length(g.cell_to_chart)

Gridap.Geometry.get_reffes(g::AtlasGrid) =
  Fill(g.reffe, Gridap.Geometry.num_cells(g))

Gridap.Geometry.get_cell_type(g::AtlasGrid) =
  Fill(Int8(1), Gridap.Geometry.num_cells(g))

Gridap.Geometry.get_cell_node_ids(g::AtlasGrid) =
  Gridap.Geometry.get_cell_node_ids(g.topo_grid)

function Gridap.Geometry.get_node_coordinates(g::AtlasGrid)
  @notimplemented """\n
  AtlasGrid is DG-style: there are no globally shared nodes.
  Use `get_cell_coordinates(g)` to access the $(num_point_dims(g))-dimensional
  corner coordinates per cell instead.
  """
end

Gridap.Geometry.get_cell_coordinates(g::AtlasGrid) = g.cell_coordinates

# ----------------------------------------------------------
# Custom API
# ----------------------------------------------------------

get_topo_grid(g::AtlasGrid)    = g.topo_grid
get_physical_maps(g::AtlasGrid) = g.physical_maps
get_cell_to_chart(g::AtlasGrid) = g.cell_to_chart
get_ambient_dim(::AtlasGrid{Dc,Da}) where {Dc,Da} = Da


# ============================================================
# AtlasDiscreteModel
# ============================================================

"""
    AtlasDiscreteModel{Dc,Da,Tp,G,A,M,O,T,L} <: Gridap.Geometry.DiscreteModel{Dc,Dc}

Combines an `AtlasGrid{Dc,Da}` (Da-dimensional geometry, DG-style) with a
`GridTopology{Dc,Dc}` (parametric topology) and a `FaceLabeling`.
"""
struct AtlasDiscreteModel{Dc, Da, Tp, G, A, M, O,
                           T <: Gridap.Geometry.GridTopology,
                           L <: Gridap.Geometry.FaceLabeling
                           } <: Gridap.Geometry.DiscreteModel{Dc,Dc}
  atlas_grid    :: AtlasGrid{Dc,Da,Tp,G,A,M,O}
  grid_topology :: T
  face_labeling :: L

  function AtlasDiscreteModel(
    atlas_grid    :: AtlasGrid{Dc,Da,Tp,G,A,M,O},
    grid_topology :: Gridap.Geometry.GridTopology{Dc,Dc},
    face_labeling :: Gridap.Geometry.FaceLabeling,
  ) where {Dc,Da,Tp,G,A,M,O}
    T = typeof(grid_topology)
    L = typeof(face_labeling)
    new{Dc,Da,Tp,G,A,M,O,T,L}(atlas_grid, grid_topology, face_labeling)
  end
end

Gridap.Geometry.get_grid(m::AtlasDiscreteModel)          = m.atlas_grid
Gridap.Geometry.get_grid_topology(m::AtlasDiscreteModel) = m.grid_topology
Gridap.Geometry.get_face_labeling(m::AtlasDiscreteModel) = m.face_labeling

get_atlas_grid(m::AtlasDiscreteModel)    = m.atlas_grid
get_ambient_dim(m::AtlasDiscreteModel{Dc,Da}) where {Dc,Da} = Da
get_physical_maps(m::AtlasDiscreteModel) = get_physical_maps(m.atlas_grid)
get_cell_to_chart(m::AtlasDiscreteModel) = get_cell_to_chart(m.atlas_grid)


# ============================================================
# Helper: transform 2D parametric → Da-dimensional ambient coordinates
# ============================================================

"""
    _transform_parametric_to_ambient(param_table, panels_local, physical_maps)

For each cell i, apply `physical_maps[panels_local[i]]` to every 2D parametric
corner in `param_table[i]`, returning a `Table{Point{Da}}` of ambient corners.

`physical_maps` must be indexable by panel id and callable on `Point{2,Float64}`.
"""
function _transform_parametric_to_ambient(
    param_table,
    panels_local::AbstractVector{Int},
    physical_maps,
)
  ncells   = length(panels_local)
  ncorners = 4   # QUAD only
  data_out = Vector{Point{3,Float64}}(undef, ncells * ncorners)
  ptrs_out = Vector{Int32}(undef, ncells + 1)
  ptrs_out[1] = 1
  for i in 1:ncells
    ptrs_out[i+1] = ptrs_out[i] + ncorners
    fwd     = physical_maps[panels_local[i]]
    corners = param_table[i]
    for j in 1:ncorners
      data_out[ptrs_out[i] + j - 1] = fwd(corners[j])
    end
  end
  Gridap.Arrays.Table(data_out, ptrs_out)
end
