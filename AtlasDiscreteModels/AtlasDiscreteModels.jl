# AtlasDiscreteModels.jl
#
# Defines AtlasDiscreteModel{Dc,Da}: wraps an AtlasGrid with topology and face labeling.
# Physical Da-dimensional coordinates are computed HERE (visualization_data only).

include("AtlasGrids.jl")   # also includes CoarseMeshes.jl

# ============================================================
# AtlasDiscreteModel
# ============================================================

"""
    AtlasDiscreteModel{Dc,Da,...} <: Gridap.Geometry.DiscreteModel{Dc,Dc}

Combines an `AtlasGrid{Dc,Da}` (local Dc-dim geometry) with a `GridTopology{Dc,Dc}`
and a `FaceLabeling`.  Physical Da-dimensional coordinates are never stored; they are
computed on the fly inside `visualization_data`.
"""
struct AtlasDiscreteModel{Dc, Da,
                           G, A, M, O,
                           T <: Gridap.Geometry.GridTopology{Dc,Dc},
                           L <: Gridap.Geometry.FaceLabeling
                           } <: Gridap.Geometry.DiscreteModel{Dc,Dc}
  atlas_grid    :: AtlasGrid{Dc,Da,G,A,M,O}
  grid_topology :: T
  face_labeling :: L

  function AtlasDiscreteModel(
    atlas_grid    :: AtlasGrid{Dc,Da,G,A,M,O},
    grid_topology :: Gridap.Geometry.GridTopology{Dc,Dc},
    face_labeling :: Gridap.Geometry.FaceLabeling,
  ) where {Dc,Da,G,A,M,O}
    T = typeof(grid_topology)
    L = typeof(face_labeling)
    new{Dc,Da,G,A,M,O,T,L}(atlas_grid, grid_topology, face_labeling)
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

get_atlas_grid(m::AtlasDiscreteModel)    = m.atlas_grid
get_ambient_dim(m::AtlasDiscreteModel{Dc,Da}) where {Dc,Da} = Da
get_physical_maps(m::AtlasDiscreteModel) = get_physical_maps(m.atlas_grid)
get_cell_to_chart(m::AtlasDiscreteModel) = get_cell_to_chart(m.atlas_grid)

# ============================================================
# Visualization: physical coords computed here only
# ============================================================

"""
    _local_to_physical(cell_local_coords, cell_to_chart, physical_maps)

Apply `physical_maps[cell_to_chart[i]]` to every local corner in `cell_local_coords[i]`,
returning a `Table` of Da-dimensional ambient corner coordinates.
Called only from `visualization_data`.

One `return_cache` is allocated per chart and reused for every corner evaluation of that
chart's map, following the Gridap cache pattern for zero-allocation inner loops.
"""
function _local_to_physical(cell_local_coords, cell_to_chart, physical_maps)
  ncells = length(cell_to_chart)

  # One cache per chart; reused across all cells/corners that share the same chart map.
  sample_pt    = cell_local_coords[1][1]
  chart_caches = [Gridap.Arrays.return_cache(physical_maps[p], sample_pt)
                  for p in eachindex(physical_maps)]
  PtOut = typeof(Gridap.Arrays.evaluate!(chart_caches[cell_to_chart[1]],
                                          physical_maps[cell_to_chart[1]], sample_pt))

  total_corners = length(cell_local_coords.data)
  data_out = Vector{PtOut}(undef, total_corners)
  ptrs_out = Vector{Int32}(undef, ncells + 1)
  ptrs_out[1] = 1

  for i in 1:ncells
    corners       = cell_local_coords[i]
    nc            = length(corners)
    ptrs_out[i+1] = ptrs_out[i] + nc
    chart         = cell_to_chart[i]
    fwd           = physical_maps[chart]
    cache         = chart_caches[chart]
    for j in eachindex(corners)
      data_out[ptrs_out[i] + j - 1] = Gridap.Arrays.evaluate!(cache, fwd, corners[j])
    end
  end

  Gridap.Arrays.Table(data_out, ptrs_out)
end

function Gridap.Visualization.visualization_data(
    model   :: AtlasDiscreteModel{Dc,Da},
    filebase :: AbstractString;
    labels  :: Gridap.Geometry.FaceLabeling = Gridap.Geometry.get_face_labeling(model),
) where {Dc,Da}
  g    = model.atlas_grid
  phys = _local_to_physical(g.cell_local_coords, g.cell_to_chart, g.physical_maps)
  node_ids = Gridap.Arrays.Table(Int32.(1:length(phys.data)), phys.ptrs)
  viz_grid = UnstructuredGrid(
    phys.data,
    node_ids,
    Gridap.Geometry.get_reffes(g.param_grid),
    Gridap.Geometry.get_cell_type(g.param_grid),
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
  coarse_model = Gridap.Geometry.UnstructuredDiscreteModel(info.grid)
  AtlasDiscreteModel(coarse_model, info.local_coords, info.physical_maps, num_refinements;
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
  coarse_model = Gridap.Geometry.UnstructuredDiscreteModel(info.grid)
  AtlasDiscreteModel(coarse_model, info.local_coords, custom_maps, num_refinements;
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
