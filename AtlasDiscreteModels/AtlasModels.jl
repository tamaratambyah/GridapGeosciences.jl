# AtlasModels.jl
#
# Serial (no-MPI) definitions of AtlasGrid, AtlasDiscreteModel, and the
# _transform_parametric_to_ambient helper.
#
# Replaces the serial portion of SBAtlasModels.jl with no parallel imports
# (GridapDistributed, GridapP4est, PartitionedArrays, MPI).

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Helpers
using Gridap.Adaptivity, Gridap.Visualization
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
# Outer constructor: build from coarse model + physical maps
# ----------------------------------------------------------

"""
    AtlasGrid(coarse_model, physical_maps, num_refinements=1; orientation_style=NonOriented())

Build an `AtlasGrid{Dc,Da}` by uniformly refining `coarse_model` `num_refinements` times
and applying `physical_maps[chart]` to each fine cell's parametric corners.

- `coarse_model`     — `DiscreteModel{Dc,Dc}` whose cells define the chart panels.
                       `physical_maps[i]` is the map for coarse cell `i`.
- `physical_maps`    — Indexable collection of callables `Point{Dc} → Point{Da}`,
                       one per coarse cell.
- `num_refinements`  — Number of uniform red-refinement passes (default 1).
"""
function AtlasGrid(
    coarse_model    :: Gridap.Geometry.DiscreteModel{Dc,Dc},
    physical_maps,
    num_refinements :: Int = 1;
    orientation_style :: Gridap.Geometry.OrientationStyle = NonOriented(),
) where Dc

  # ── Iterative refinement, accumulating fine→coarse-panel mapping ──────────
  current_model = coarse_model
  # cell_to_chart[i] = coarse panel index for cell i of current_model
  cell_to_chart = collect(1:num_cells(coarse_model))

  for _ in 1:num_refinements
    adapted       = Gridap.Adaptivity.refine(current_model)
    level_map     = adapted.glue.n2o_faces_map[Dc+1]   # fine cell → previous cell
    cell_to_chart = cell_to_chart[level_map]            # compose: fine → original coarse
    current_model = adapted.model
  end

  fine_model = current_model
  fine_grid  = Gridap.Geometry.get_grid(fine_model)

  # ── Reference element from the fine grid ─────────────────────────────────
  reffes   = Gridap.Geometry.get_reffes(fine_grid)
  reffe    = reffes[1]   # all cells share the same reffe (uniform mesh)

  # ── Parametric corner coordinates of fine cells ───────────────────────────
  param_coords = Gridap.Geometry.get_cell_coordinates(fine_grid)

  # ── Map parametric corners to ambient space ───────────────────────────────
  cell_coords_da = _transform_parametric_to_ambient(
    param_coords, cell_to_chart, physical_maps)

  AtlasGrid(fine_grid, cell_coords_da, cell_to_chart, physical_maps, reffe, orientation_style)
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

# DG-sequential node IDs: cell i owns nodes ptrs[i]..ptrs[i+1]-1 in the flat coord array.
function Gridap.Geometry.get_cell_node_ids(g::AtlasGrid)
  cc = g.cell_coordinates   # Table{Point{Da}}; cc.ptrs gives per-cell offsets
  Gridap.Arrays.Table(Int32.(1:length(cc.data)), cc.ptrs)
end

# The flat DG coordinate array: all cell corners laid out sequentially.
Gridap.Geometry.get_node_coordinates(g::AtlasGrid) = g.cell_coordinates.data

Gridap.Geometry.get_cell_coordinates(g::AtlasGrid) = g.cell_coordinates

# ----------------------------------------------------------
# Custom API
# ----------------------------------------------------------

get_topo_grid(g::AtlasGrid)     = g.topo_grid
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

# ----------------------------------------------------------
# Outer constructor: build from coarse model + physical maps
# ----------------------------------------------------------

"""
    AtlasDiscreteModel(coarse_model, physical_maps, num_refinements=1; orientation_style=NonOriented())

Convenience constructor: refines `coarse_model`, builds an `AtlasGrid`, then
wraps it with the topology and labeling from the finest refined model.
"""
function AtlasDiscreteModel(
    coarse_model    :: Gridap.Geometry.DiscreteModel{Dc,Dc},
    physical_maps,
    num_refinements :: Int = 1;
    orientation_style :: Gridap.Geometry.OrientationStyle = NonOriented(),
) where Dc

  # Refine to get topology and labeling of the fine model
  current_model = coarse_model
  for _ in 1:num_refinements
    current_model = Gridap.Adaptivity.refine(current_model).model
  end
  fine_model = current_model

  atlas_grid = AtlasGrid(coarse_model, physical_maps, num_refinements;
                         orientation_style)

  AtlasDiscreteModel(
    atlas_grid,
    Gridap.Geometry.get_grid_topology(fine_model),
    Gridap.Geometry.get_face_labeling(fine_model),
  )
end

Gridap.Geometry.get_grid(m::AtlasDiscreteModel)          = m.atlas_grid
Gridap.Geometry.get_grid_topology(m::AtlasDiscreteModel) = m.grid_topology
Gridap.Geometry.get_face_labeling(m::AtlasDiscreteModel) = m.face_labeling

get_atlas_grid(m::AtlasDiscreteModel)    = m.atlas_grid
get_ambient_dim(m::AtlasDiscreteModel{Dc,Da}) where {Dc,Da} = Da
get_physical_maps(m::AtlasDiscreteModel) = get_physical_maps(m.atlas_grid)
get_cell_to_chart(m::AtlasDiscreteModel) = get_cell_to_chart(m.atlas_grid)


# ============================================================
# Helper: transform parametric corners → ambient coordinates
# ============================================================

"""
    _transform_parametric_to_ambient(param_coords, cell_to_chart, physical_maps)

For each cell `i`, apply `physical_maps[cell_to_chart[i]]` to every parametric
corner in `param_coords[i]`, returning a `Table` of ambient-space corners.

The number of corners per cell is read from `param_coords[i]` directly, so this
works for any polytope (QUAD, HEX, TRI, …). The ambient `Point` type is inferred
from the first map evaluation, so no dimension is hardcoded.
"""
function _transform_parametric_to_ambient(
    param_coords,
    cell_to_chart :: AbstractVector{Int},
    physical_maps,
)
  ncells = length(cell_to_chart)

  # Infer output Point type from first map evaluation
  sample_pt  = param_coords[1][1]
  sample_out = physical_maps[cell_to_chart[1]](sample_pt)
  PtOut      = typeof(sample_out)

  # Total number of corners across all cells (variable per cell in general)
  total_corners = sum(length(param_coords[i]) for i in 1:ncells)

  data_out = Vector{PtOut}(undef, total_corners)
  ptrs_out = Vector{Int32}(undef, ncells + 1)
  ptrs_out[1] = 1

  for i in 1:ncells
    corners  = param_coords[i]
    nc       = length(corners)
    ptrs_out[i+1] = ptrs_out[i] + nc
    fwd      = physical_maps[cell_to_chart[i]]
    for j in eachindex(corners)
      data_out[ptrs_out[i] + j - 1] = fwd(corners[j])
    end
  end

  Gridap.Arrays.Table(data_out, ptrs_out)
end


# ============================================================
# Visualization (writevtk support)
# ============================================================

# AtlasGrid now implements get_node_coordinates / get_cell_node_ids with DG-sequential
# semantics, so write_vtk_file(atlas_grid, ...) works without any extra grid object.
#
# For AtlasDiscreteModel we only write the cell-dimension (d=Dc) file; the
# sub-dimensional face grids (edges, vertices) would need parametric coordinates
# from grid_topology, which is a separate concern.

function Gridap.Visualization.visualization_data(
    model::AtlasDiscreteModel{Dc,Da}, filebase::AbstractString;
    labels::Gridap.Geometry.FaceLabeling = Gridap.Geometry.get_face_labeling(model),
) where {Dc,Da}
  Gridap.Visualization.visualization_data(model.atlas_grid, "$(filebase)_$(Dc)")
end
