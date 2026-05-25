using FillArrays
using Gridap
using GridapGeosciences

using MPI
using PartitionedArrays
using GridapDistributed
using GridapP4est
using P4est_wrapper

include("AtlasGrid.jl")

"""
Bilinear map from p4est integer quadrant coordinates (x, y) to physical space,
given the four corner vertices of the coarse tree in p4est ordering:
  v0: (xi=0, yi=0),  v1: (xi=1, yi=0)
  v2: (xi=0, yi=1),  v3: (xi=1, yi=1)

Equivalent to `p4est_qcoord_to_vertex` (2D, P4_TO_P8=false) but takes
tree vertices directly instead of a connectivity pointer + treeid.
"""
function quad_coord_to_vertex_coord(coarse_vertices,
                                    x::p4est_qcoord_t,
                                    y::p4est_qcoord_t)
  P4EST_ROOT_LEN = p4est_qcoord_t(1) << P4est_wrapper.P4EST_MAXLEVEL

  tx = Float64(x) / Float64(P4EST_ROOT_LEN)
  ty = Float64(y) / Float64(P4EST_ROOT_LEN)

  wx = (1.0 - tx, tx)
  wy = (1.0 - ty, ty)

  vx = 0.0
  vy = 0.0
  k = 1
  for yi in 1:2
    for xi in 1:2
      w = wy[yi] * wx[xi]
      vx += w * coarse_vertices[k][1]
      vy += w * coarse_vertices[k][2]
      k += 1
    end
  end
  vx, vy
end

function get_quad_vertex_coord(coarse_vertices,
                               x::p4est_qcoord_t,
                               y::p4est_qcoord_t,
                               level::Int8,
                               corner::Cint)

  myself = Ref{p4est_quadrant_t}(
    p4est_quadrant_t(x,y,level,Int8(0),Int16(0),P4est_wrapper.quadrant_data(Clong(0)))
  )
  neighbour = Ref{p4est_quadrant_t}(myself[])
  if corner == 1
      p4est_quadrant_face_neighbor(myself,corner,neighbour)
  elseif corner == 2
      p4est_quadrant_face_neighbor(myself,corner+1,neighbour)
  elseif corner == 3
      p4est_quadrant_corner_neighbor(myself,corner,neighbour)
  end
  # Extract numerical coordinates of lower_left
  # corner of my corner neighbour
  quad_coord_to_vertex_coord(coarse_vertices,
                             neighbour[].x,
                             neighbour[].y)
end


function setup_atlas_coarse_model_original()
    #        x=-2                                    x=2
    # y=2     5───────────────────────────────────────8
    #         │╲                                     ╱│
    #         │  ╲              C4                /   │
    # y=1     │    1───────────────────────────4      │
    #         │    │                           │      │
    #         │ C5 │y          C1              │ C6   │
    #         │    │. x                        │      │
    # y=-1    │    2───────────────────────────3      │
    #         │. ╱              C3               ╲    │
    #         │╱.                                  ╲ .│
    # y=-2    6───────────────────────────────────────7

    # C1: [2,3,1,4] — inner square
    # C2: [6,7,5,8] — outer square
    # C3: [6,7,2,3] — bottom trapezoid
    # C4: [5,1,8,4] — top trapezoid
    # C5: [6,2,5,1] — left trapezoid
    # C6: [7,8,3,4] — right trapezoid

    ptr  = [ 1, 5, 9, 13, 17, 21, 25 ]
    data = [ 2,3,1,4, 6,7,5,8, 6,7,2,3, 5,1,8,4, 6,2,5,1, 7,8,3,4 ]
    
    cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
    # Atlas coordinates. They are not actually needed.
    node_coordinates = Vector{Point{2,Float64}}(undef,8)
    node_coordinates[1]=Point{2,Float64}(-1.0,  1.0)
    node_coordinates[2]=Point{2,Float64}(-1.0, -1.0)
    node_coordinates[3]=Point{2,Float64}( 1.0, -1.0)
    node_coordinates[4]=Point{2,Float64}( 1.0,  1.0)
    node_coordinates[5]=Point{2,Float64}(-2.0,  2.0)
    node_coordinates[6]=Point{2,Float64}(-2.0, -2.0)
    node_coordinates[7]=Point{2,Float64}( 2.0, -2.0)
    node_coordinates[8]=Point{2,Float64}( 2.0,  2.0)

    polytope=QUAD
    scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
    cell_types=collect(Fill(1,length(cell_vertex_lids)))
    cell_reffes=[scalar_reffe]
    grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                            cell_vertex_lids,
                                            cell_reffes,
                                            cell_types,
                                            Gridap.Geometry.Oriented())
    Gridap.Geometry.UnstructuredDiscreteModel(grid)
  end

  function setup_atlas_coarse_model_santi_version_of_tamara()
    #        x=-2                                    x=2
    # y=2     5───────────────────────────────────────8
    #         │╲                                     ╱│
    #         │  ╲              C4                /   │
    # y=1     │    1───────────────────────────4      │
    #         │    │                           │      │
    #         │ C5 │y          C1              │ C6   │
    #         │    │. x                        │      │
    # y=-1    │    2───────────────────────────3      │
    #         │. ╱              C3               ╲    │
    #         │╱.                                  ╲ .│
    # y=-2    6───────────────────────────────────────7

    # C1: [2,3,1,4] — inner square
    # C2: [6,7,5,8] — outer square
    # C3: [6,7,2,3] — bottom trapezoid
    # C4: [5,1,8,4] — top trapezoid
    # C5: [6,2,5,1] — left trapezoid
    # C6: [7,8,3,4] — right trapezoid

    ptr  = [ 1, 5, 9, 13, 17, 21, 25 ]
    data = [ 2,3,1,4, 6,7,5,8, 6,7,2,3, 5,1,8,4, 6,2,5,1, 7,8,3,4 ]
    
    cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
    # Atlas coordinates. They are not actually needed.
    node_coordinates = Vector{Point{2,Float64}}(undef,8)
    node_coordinates[1]=Point{2,Float64}(-1.0, -1.0)
    node_coordinates[2]=Point{2,Float64}( 1.0, -1.0)
    node_coordinates[3]=Point{2,Float64}(-2.0, -2.0)
    node_coordinates[4]=Point{2,Float64}( 2.0, -2.0)
    node_coordinates[5]=Point{2,Float64}(-2.0,  2.0)
    node_coordinates[6]=Point{2,Float64}( 2.0,  2.0)
    node_coordinates[7]=Point{2,Float64}( 1.0,  1.0)
    node_coordinates[8]=Point{2,Float64}(-1.0,  1.0)

    polytope=QUAD
    scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
    cell_types=collect(Fill(1,length(cell_vertex_lids)))
    cell_reffes=[scalar_reffe]
    grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                            cell_vertex_lids,
                                            cell_reffes,
                                            cell_types,
                                            Gridap.Geometry.Oriented())
    Gridap.Geometry.UnstructuredDiscreteModel(grid)
  end

  function setup_atlas_coarse_model_alberto_version_of_tamara()
    #        x=-2                                    x=2
    # y=2     1───────────────────────────────────────8  C5
    #         │╲                                     ╱│
    #         │  ╲              C6                /   │
    # y=1     │    3───────────────────────────5      │
    #         │    │.y                         │      │
    #         │ C1 │x           C2             │ C4   │
    #         │    │                           │      │
    # y=-1    │    4───────────────────────────6      │
    #         │  ╱              C3               ╲    │
    #         │╱                                   ╲  │
    # y=-2    2───────────────────────────────────────7

    # C1: [1,2,3,4] — left trapezoid
    # C2: [3,4,5,6] — inner square
    # C3: [2,7,4,6] — bottom trapezoid
    # C4: [8,5,7,6] — right trapezoid
    # C5: [1,8,2,7] — outer square
    # C6: [1,3,8,5] — top trapezoid

    ptr  = [ 1, 5, 9, 13, 17, 21, 25 ]
    data = [ 1,2,3,4, 3,4,5,6,  2,7,4,6, 8,5,7,6, 1,8,2,7, 1,3,8,5 ]
    
    cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
    # Atlas coordinates. They are not actually needed.
    node_coordinates = Vector{Point{2,Float64}}(undef,8)
    node_coordinates[1]=Point{2,Float64}(-2.0,  2.0)
    node_coordinates[2]=Point{2,Float64}(-2.0, -2.0)
    node_coordinates[3]=Point{2,Float64}(-1.0,  1.0)
    node_coordinates[4]=Point{2,Float64}(-1.0, -1.0)
    node_coordinates[5]=Point{2,Float64}( 1.0,  1.0)
    node_coordinates[6]=Point{2,Float64}( 1.0, -1.0)
    node_coordinates[7]=Point{2,Float64}( 2.0, -2.0)
    node_coordinates[8]=Point{2,Float64}( 2.0,  2.0)

    polytope=QUAD
    scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
    cell_types=collect(Fill(1,length(cell_vertex_lids)))
    cell_reffes=[scalar_reffe]
    grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                            cell_vertex_lids,
                                            cell_reffes,
                                            cell_types,
                                            Gridap.Geometry.Oriented())
    Gridap.Geometry.UnstructuredDiscreteModel(grid)
  end


# DEPRECATED, not used anymore
# function setup_chart_maps_to_cube_surface(coarse_model::DiscreteModel)
#     polytope=QUAD
#     scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
#     shape_funs = Gridap.ReferenceFEs.get_shapefuns(scalar_reffe)
#     cell_shape_funs = Fill(shape_funs, num_cells(coarse_model))

#     # 3D cube corner nodes shared by all six chart maps
#     cube_nodes = Vector{Point{3,Float64}}(undef,8)
#     cube_nodes[1]=Point{3,Float64}(-1.0, -1.0,  1.0)
#     cube_nodes[2]=Point{3,Float64}( 1.0, -1.0,  1.0)
#     cube_nodes[3]=Point{3,Float64}( 1.0,  1.0,  1.0)
#     cube_nodes[4]=Point{3,Float64}(-1.0,  1.0,  1.0)
#     cube_nodes[5]=Point{3,Float64}(-1.0, -1.0, -1.0)
#     cube_nodes[6]=Point{3,Float64}( 1.0, -1.0, -1.0)
#     cube_nodes[7]=Point{3,Float64}( 1.0,  1.0, -1.0)
#     cube_nodes[8]=Point{3,Float64}(-1.0,  1.0, -1.0)

#     cell_node_ids = Gridap.Geometry.get_cell_node_ids(coarse_model)

#     # Maps from reference [0,1]² → 3D cube face (bilinear, one per coarse cell)
#     m3d = Broadcasting(Reindex(cube_nodes))
#     cell_cube_node_coords = lazy_map(m3d, cell_node_ids)
#     ref_to_cube_maps = lazy_map(Gridap.Fields.linear_combination, cell_cube_node_coords, cell_shape_funs)
    
#     # Maps from reference [0,1]² → 2D physical atlas domain 
#     # (from coarse grid geometry)
#     ref_to_atlas_maps = Gridap.Geometry.get_cell_map(get_grid(coarse_model))

#     # Compose: 2D physical atlas → 3D cube face
#     # = (ref→cube) ∘ inverse(ref→atlas)
#     atlas_to_ref_maps = lazy_map(Gridap.Fields.inverse_map, ref_to_atlas_maps)
#     lazy_map(∘, ref_to_cube_maps, atlas_to_ref_maps)
# end

# DEPRECATED, not used anymore
# function generate_cube_surface_dmodel(octree_dmodel::OctreeDistributedDiscreteModel{2,2}, cell_to_cube_surface_map)
#     ptr_pXest_ghost=GridapP4est.setup_pXest_ghost(octree_dmodel.pXest_type,
#                                                   octree_dmodel.ptr_pXest)
    
#     cell_coordinates = 
#        generate_cell_coordinates(octree_dmodel.parts,
#                              [ CUBE_HALF_EDGE.*[Point(-1.0,-1.0),Point(1.0,-1.0),Point(-1.0,1.0),Point(1.0,1.0)]
#                              for i=1:6],
#                               octree_dmodel.ptr_pXest_connectivity,
#                               octree_dmodel.ptr_pXest,
#                               ptr_pXest_ghost)
#    lmodels = map(local_views(octree_dmodel.dmodel), 
#                  cell_to_cube_surface_map,
#                  cell_coordinates) do omodel, cell_to_cube_surface_map, cell_coordinates
#        ogrid = get_grid(omodel)
#        cell_node_coords = get_cell_coordinates(ogrid)
#        cell_to_nodes = Gridap.Geometry.get_cell_node_ids(ogrid)
                                                   
#        surface_cell_node_coordinates = lazy_map(evaluate, cell_to_cube_surface_map, cell_coordinates)

#        eT = eltype(Gridap.Arrays.testitem(surface_cell_node_coordinates))
#        surface_nodes = zeros(eT, length(Gridap.Geometry.get_node_coordinates(ogrid)))
#        Gridap.FESpaces._free_and_dirichlet_values_fill!(
#                 surface_nodes, eT[],
#                 array_cache(surface_cell_node_coordinates),
#                 array_cache(cell_to_nodes),
#                 surface_cell_node_coordinates,
#                 cell_to_nodes,
#                 eachindex(cell_to_nodes)
#        ) 
#        surface_grid = Gridap.Geometry.UnstructuredGrid(surface_nodes,
#                         cell_to_nodes,
#                         Gridap.Geometry.get_reffes(ogrid),
#                         Gridap.Geometry.get_cell_type(ogrid),
#                         Gridap.Geometry.OrientationStyle(ogrid))

#         topo = Gridap.Geometry.UnstructuredGridTopology(surface_nodes,
#                                                         cell_to_nodes,
#                                                         Gridap.Geometry.get_cell_type(ogrid),
#                                                         Gridap.Geometry.get_polytopes(ogrid),
#                                                         Gridap.Geometry.OrientationStyle(ogrid))
#         labels = get_face_labeling(omodel)
#         Gridap.Geometry.UnstructuredDiscreteModel(surface_grid, topo, labels)                
#    end   
#    GridapDistributed.GenericDistributedDiscreteModel(lmodels, GridapDistributed.get_cell_gids(octree_dmodel.dmodel))
# end


function generate_cell_coordinates(ranks,
                                   coarse_cell_wise_vertex_coordinates,
                                   ptr_pXest_connectivity,
                                   ptr_pXest,
                                   ptr_pXest_ghost)
  Dc=2
  PXEST_CORNERS=4
  pXest_ghost = ptr_pXest_ghost[]
  pXest = ptr_pXest[]

  # Obtain ghost quadrants
  ptr_ghost_quadrants = Ptr{P4est_wrapper.p4est_quadrant_t}(pXest_ghost.ghosts.array)

  tree_offsets = unsafe_wrap(Array, pXest_ghost.tree_offsets, pXest_ghost.num_trees+1)
  map(ranks) do _
     ncells=pXest.local_num_quadrants+pXest_ghost.ghosts.elem_count
     panels = Vector{Int}(undef,ncells)
     data = Vector{Point{Dc,Float64}}(undef,ncells*PXEST_CORNERS)
     ptr  = GridapGeosciences.generate_ptr(Dc,ncells)
     current=1
     for itree=1:pXest_ghost.num_trees
       tree = p4est_tree_array_index(pXest.trees, itree-1)[]
       for cell=1:tree.quadrants.elem_count
          quadrant=p4est_quadrant_array_index(tree.quadrants, cell-1)[]
          for vertex=1:PXEST_CORNERS
            vxy = get_quad_vertex_coord(coarse_cell_wise_vertex_coordinates[itree],
                                        quadrant.x,
                                        quadrant.y,
                                        quadrant.level,
                                        Cint(vertex-1))
            data[current]=Point{Dc,Float64}(vxy...)
            current=current+1
          end
       end
     end

     # Go over ghost cells
     for i=1:pXest_ghost.num_trees
       for j=tree_offsets[i]:tree_offsets[i+1]-1
          panels[current_cell]=i
          quadrant = ptr_ghost_quadrants[j+1]
          for vertex=1:PXEST_CORNERS
            vxy = get_quad_vertex_coord(coarse_cell_wise_vertex_coordinates[i],
                                        quadrant.x,
                                        quadrant.y,
                                        quadrant.level,
                                        Cint(vertex-1))
            data[current]=Point{Dc,Float64}(vxy...)
            current=current+1
         end
       end
     end
     Gridap.Arrays.Table(data,ptr)
  end
end

function generate_cell_to_chart_id(ranks,
                                   coarse_discrete_model,
                                   ptr_pXest_connectivity,
                                   ptr_pXest,
                                   ptr_pXest_ghost)

  pXest_ghost = ptr_pXest_ghost[]
  pXest = ptr_pXest[]

  coarse_cell_to_chart_id = collect(1:num_cells(coarse_discrete_model))

  tree_offsets = unsafe_wrap(Array, pXest_ghost.tree_offsets, pXest_ghost.num_trees+1)
  cell_to_chart_id=map(ranks) do part
     ncells=pXest.local_num_quadrants+pXest_ghost.ghosts.elem_count
     cell_to_chart_id = Vector{Int}(undef,ncells)
     current_cell=1
     for itree=1:pXest_ghost.num_trees
       tree = GridapP4est.p4est_tree_array_index(pXest.trees, itree-1)[]
       for cell=1:tree.quadrants.elem_count
          cell_to_chart_id[current_cell]=coarse_cell_to_chart_id[itree]
          current_cell=current_cell+1
       end
     end
     # Go over ghost cells
     for i=1:pXest_ghost.num_trees
      for j=tree_offsets[i]:tree_offsets[i+1]-1
          cell_to_chart_id[current_cell]=i
          current_cell=current_cell+1
       end
     end
     cell_to_chart_id
  end
  cell_to_chart_id
end

function generate_cell_to_chart_id(octree_dmodel::OctreeDistributedDiscreteModel)
   pXest_type = octree_dmodel.pXest_type
   ptr_pXest = octree_dmodel.ptr_pXest
   # Build the ghost layer
   ptr_pXest_ghost=GridapP4est.setup_pXest_ghost(pXest_type,ptr_pXest)
   cell_to_chart_id = 
      generate_cell_to_chart_id(octree_dmodel.parts,
                                octree_dmodel.coarse_model,
                                octree_dmodel.ptr_pXest_connectivity,
                                octree_dmodel.ptr_pXest,
                                ptr_pXest_ghost)
   GridapP4est.pXest_ghost_destroy(pXest_type,ptr_pXest_ghost)
   cell_to_chart_id
end

function generate_atlas_grids(octree_dmodel::OctreeDistributedDiscreteModel)
   pXest_ghost = GridapP4est.setup_pXest_ghost(octree_dmodel.pXest_type,octree_dmodel.ptr_pXest)
   cell_coordinates = 
      generate_cell_coordinates(octree_dmodel.parts,
                                [GridapGeosciences.Geometry.CUBE_HALF_EDGE.*[Point(-1.0,-1.0),Point(1.0,-1.0),Point(-1.0,1.0),Point(1.0,1.0)] for i=1:6 ],
                                octree_dmodel.ptr_pXest_connectivity,
                                octree_dmodel.ptr_pXest,
                                pXest_ghost
                               )
   GridapP4est.pXest_ghost_destroy(octree_dmodel.pXest_type,pXest_ghost)
   cell_to_chart_id = generate_cell_to_chart_id(octree_dmodel)
   reffe = Gridap.ReferenceFEs.ReferenceFE(QUAD,Gridap.ReferenceFEs.lagrangian,Float64,1)
   orientation_style = Gridap.Geometry.Oriented()
   map(cell_coordinates, cell_to_chart_id, local_views(octree_dmodel.dmodel)) do cell_coordinates, 
                                                                                 cell_to_chart_id, 
                                                                                 model
      AtlasGrid(cell_coordinates,
                Gridap.Geometry.get_cell_node_ids(get_grid(model)),
                cell_to_chart_id,
                reffe,
                orientation_style)
   end
end 

struct AtlasOctreeDistributedDiscreteModel{A<:OctreeDistributedDiscreteModel{2,2},B} <: GridapDistributed.DistributedDiscreteModel{2,2}
  octree_dmodel::A
  atlas_grids::B
end

# cell_to_chart_id = generate_cell_to_chart_id(octree_dmodel)

function AtlasOctreeDistributedDiscreteModel(ranks, num_uniform_refinements)
  coarse_model = setup_atlas_coarse_model_alberto_version_of_tamara()
  octree_dmodel = OctreeDistributedDiscreteModel(ranks, coarse_model, num_uniform_refinements)
  atlas_grids = generate_atlas_grids(octree_dmodel)
  AtlasOctreeDistributedDiscreteModel(octree_dmodel, atlas_grids)
end 

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
ranks = distribute_with_mpi(LinearIndices((prod(nprocs),)))
model = AtlasOctreeDistributedDiscreteModel(ranks, 1)

model_old = CubedSphere2DParametricOctreeDistributedDiscreteModel(ranks,1.0;num_initial_uniform_refinements=1)

cold = model_old.parametric_dmodel.models.item_ref[].grid.cell_map.args[1]
cnew = model.atlas_grids.item_ref[].cell_coordinates
for (co,cn) in zip(cold, cnew)
   @assert co ≈ cn
end


