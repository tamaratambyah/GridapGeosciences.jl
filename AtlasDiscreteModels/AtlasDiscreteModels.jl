# AtlasDiscreteModels.jl
#
# Defines AtlasDiscreteModel{Dc,Da}: wraps an AtlasGrid with topology and face labeling.
# Physical Da-dimensional coordinates are computed HERE (visualization_data only).

include("AtlasGrids.jl")   # also includes CoarseMeshes.jl

# ============================================================
# AtlasDiscreteModel
# ============================================================

"""
    AtlasDiscreteModel{Dc,Da,G,A,P,O,T,L} <: Gridap.Geometry.DiscreteModel{Dc,Dc}

Combines an `AtlasGrid{Dc,Da}` (local Dc-dim geometry) with a `GridTopology{Dc,Dc}`
and a `FaceLabeling`.  Physical Da-dimensional coordinates are never stored; they are
computed on the fly inside `visualization_data`.

Face labels from the coarse `CoarseMeshInfo` model (e.g. "bottom", "top" for the
cylinder) are propagated to the fine mesh by Gridap's refinement machinery and are
accessible via `Gridap.Geometry.get_face_labeling(model)`.
"""
struct AtlasDiscreteModel{Dc, Da,
                           G, A, P, O,
                           T <: Gridap.Geometry.GridTopology{Dc,Dc},
                           L <: Gridap.Geometry.FaceLabeling
                           } <: Gridap.Geometry.DiscreteModel{Dc,Dc}
  atlas_grid    :: AtlasGrid{Dc,Da,G,A,P,O}
  grid_topology :: T
  face_labeling :: L

  function AtlasDiscreteModel(
    atlas_grid    :: AtlasGrid{Dc,Da,G,A,P,O},
    grid_topology :: Gridap.Geometry.GridTopology{Dc,Dc},
    face_labeling :: Gridap.Geometry.FaceLabeling,
  ) where {Dc,Da,G,A,P,O}
    T = typeof(grid_topology)
    L = typeof(face_labeling)
    new{Dc,Da,G,A,P,O,T,L}(atlas_grid, grid_topology, face_labeling)
  end
end

# ----------------------------------------------------------
# Outer constructor
# ----------------------------------------------------------

"""
    AtlasDiscreteModel(coarse_model, coarse_local_coords, physical_maps, num_refinements=1;
                       orientation_style=nothing)

Refine `coarse_model` `num_refinements` times, build an `AtlasGrid` with local coords,
and wrap it with the fine model's topology and face labeling.
"""
function AtlasDiscreteModel(
    coarse_model    :: Gridap.Geometry.DiscreteModel{Dc,Dc},
    coarse_local_coords,
    physical_maps,
    num_refinements :: Int = 1;
    orientation_style = nothing,
) where Dc
  # _build_atlas_grid refines exactly once and returns both the atlas_grid and the
  # fine_model, so topology and face labeling are extracted without a second refinement.
  atlas_grid, fine_model = _build_atlas_grid(
    coarse_model, coarse_local_coords, physical_maps, num_refinements, orientation_style)

  AtlasDiscreteModel(
    atlas_grid,
    Gridap.Geometry.get_grid_topology(fine_model),
    Gridap.Geometry.get_face_labeling(fine_model),
  )
end

# ----------------------------------------------------------
# DiscreteModel{Dc,Dc} interface
# ----------------------------------------------------------

Gridap.Geometry.get_grid(m::AtlasDiscreteModel)          = m.atlas_grid
Gridap.Geometry.get_grid_topology(m::AtlasDiscreteModel) = m.grid_topology
Gridap.Geometry.get_face_labeling(m::AtlasDiscreteModel) = m.face_labeling

# ----------------------------------------------------------
# Custom API
# ----------------------------------------------------------

get_atlas_grid(m::AtlasDiscreteModel)              = m.atlas_grid
get_ambient_dim(m::AtlasDiscreteModel{Dc,Da}) where {Dc,Da} = Da
get_cell_physical_maps(m::AtlasDiscreteModel)      = get_cell_physical_maps(m.atlas_grid)

# ============================================================
# Visualization: physical coords computed here only
# ============================================================

"""
    _local_to_physical(cell_local_coords, cell_physical_maps)

Return a lazy array whose `i`-th entry is the `Da`-dimensional physical corners of
cell `i`, obtained by applying `cell_physical_maps[i]` pointwise to
`cell_local_coords[i]`.

`cell_physical_maps[i]` is a `Point{Dc} → Point{Da}` map (e.g. `ForwardMap2D`).
`Broadcasting(cell_physical_maps[i])` lifts it to `Vector{Point} → Vector{Point}`;
`lazy_map(evaluate, …)` chains the two lazy arrays with zero allocation until accessed.
"""
function _local_to_physical(cell_local_coords, cell_physical_maps)
  cell_maps = lazy_map(Broadcasting, cell_physical_maps)
  lazy_map(evaluate, cell_maps, cell_local_coords)
end

function Gridap.Visualization.visualization_data(
    model   :: AtlasDiscreteModel{Dc,Da},
    filebase :: AbstractString;
    labels  :: Gridap.Geometry.FaceLabeling = Gridap.Geometry.get_face_labeling(model),
) where {Dc,Da}
  g         = model.atlas_grid
  phys_lazy = _local_to_physical(g.cell_local_coords, g.cell_physical_maps)
  ncells    = Gridap.Geometry.num_cells(g)
  n_corners = length(g.cell_local_coords[1])
  dg_node_ids = Gridap.Arrays.Table(
    Int32.(1:ncells*n_corners),
    Int32[1 + i*n_corners for i in 0:ncells],
  )
  viz_grid  = UnstructuredGrid(
    collect(Iterators.flatten(phys_lazy)),
    dg_node_ids,
    get_reffes(g),
    get_cell_type(g),
    NonOriented(),
  )
  Gridap.Visualization.visualization_data(viz_grid, "$(filebase)_$(Dc)")
end

# ============================================================
# Convenience constructors via CoarseMeshInfo / AbstractCoarseMesh
# ============================================================

# AtlasGrid equivalents live in AtlasGrids.jl (after include("CoarseMeshes.jl")).

# ----------------------------------------------------------
# AtlasDiscreteModel
# ----------------------------------------------------------

"""
    AtlasDiscreteModel(info::CoarseMeshInfo, num_refinements; orientation_style=nothing)

Build an `AtlasDiscreteModel` from a `CoarseMeshInfo`, using the physical maps stored in `info`.
"""
function AtlasDiscreteModel(
    info            :: CoarseMeshInfo,
    num_refinements :: Int;
    orientation_style = nothing,
)
  AtlasDiscreteModel(info.model, info.local_coords, info.physical_maps, num_refinements;
                     orientation_style)
end

"""
    AtlasDiscreteModel(info::CoarseMeshInfo, custom_maps, num_refinements; orientation_style=nothing)

Build an `AtlasDiscreteModel` from a `CoarseMeshInfo`, overriding the default physical maps.
"""
function AtlasDiscreteModel(
    info            :: CoarseMeshInfo,
    custom_maps,
    num_refinements :: Int;
    orientation_style = nothing,
)
  AtlasDiscreteModel(info.model, info.local_coords, custom_maps, num_refinements;
                     orientation_style)
end

"""
    AtlasDiscreteModel(mesh::AbstractCoarseMesh, num_refinements; orientation_style=nothing)

Build an `AtlasDiscreteModel` directly from a mesh descriptor (e.g. `CubedSphereMesh(1.0)`).
Calls `get_coarse_mesh(mesh)` internally and uses the default physical maps.
"""
function AtlasDiscreteModel(
    mesh            :: AbstractCoarseMesh,
    num_refinements :: Int;
    orientation_style = nothing,
)
  AtlasDiscreteModel(get_coarse_mesh(mesh), num_refinements; orientation_style)
end
