using Gridap
include("1_cell_per_panel.jl")
include("3_cells.jl")

cell_vertices_coordinates = get_cell_coordinates(cube_grid)

cell_vertices=Gridap.Geometry.get_faces(topo,2,0)
cell_edges=Gridap.Geometry.get_faces(topo,2,1)
edge_cells=Gridap.Geometry.get_faces(topo,1,2)
edge_vertices=Gridap.Geometry.get_faces(topo,1,0)

edge_vertices_sorted=deepcopy(collect(edge_vertices))
for edge=1:length(edge_vertices_sorted)
    sort!(edge_vertices_sorted[edge])
end


function get_next_unoriented_cell(edge_orientation,cell_edges,current_cell)
    for cell=current_cell:length(cell_edges)
    for e=1:length(cell_edges[cell])
          edge=cell_edges[cell][e]
          if !haskey(edge_orientation,edge)
            return cell
          end
    end
    end
    return -1
end

function orient_one_set_of_parallel_edges!(edge_orientation,
                                          cell_vertices,
                                          cell_edges,
                                          edge_cells,
                                          edge_vertices,
                                          current_cell,
                                          lid_starting_edge)
    quad_edge_vertices=Gridap.Geometry.get_face_vertices(QUAD,1)
    gid_starting_edge=cell_edges[current_cell][lid_starting_edge]
    if (edge_vertices[gid_starting_edge][1]==cell_vertices[current_cell][quad_edge_vertices[lid_starting_edge][1]])
      @assert edge_vertices[gid_starting_edge][2]==cell_vertices[current_cell][quad_edge_vertices[lid_starting_edge][2]]
      edge_orientation[gid_starting_edge]=BACKWARD
    else
      @assert edge_vertices[gid_starting_edge][1]==cell_vertices[current_cell][quad_edge_vertices[lid_starting_edge][2]]
      @assert edge_vertices[gid_starting_edge][2]==cell_vertices[current_cell][quad_edge_vertices[lid_starting_edge][1]]
      edge_orientation[gid_starting_edge]=FORWARD
    end

    δk=Set{Int}()
    δk_minus_one=Set(gid_starting_edge)
    while (length(δk_minus_one)>0)
      empty!(δk)
      for gid_edge in δk_minus_one
          for cell_edge in edge_cells[gid_edge]
            lid_cell_edge=findfirst(l->(l==gid_edge),cell_edges[cell_edge])
            first_edge_vertex=edge_orientation[gid_edge]==FORWARD ?
                                                edge_vertices[gid_edge][1] :
                                                edge_vertices[gid_edge][2]
            first_edge_vertex_cell_edge=cell_vertices[cell_edge][quad_edge_vertices[lid_cell_edge][1]]
            @assert first_edge_vertex==first_edge_vertex_cell_edge ||
                      first_edge_vertex==cell_vertices[cell_edge][quad_edge_vertices[lid_cell_edge][2]]


            opposite_edge = cell_edges[cell_edge][parallel_edges[lid_cell_edge]]
            if first_edge_vertex==first_edge_vertex_cell_edge
                o=1
            else
                o=2
            end
            first_opposite_edge_vertex=
            cell_vertices[cell_edge][quad_edge_vertices[parallel_edges[lid_cell_edge]][o]]

            if edge_vertices[opposite_edge][1]==first_opposite_edge_vertex
                opposite_edge_orientation=FORWARD
            else
                opposite_edge_orientation=BACKWARD
            end

            if !haskey(edge_orientation,opposite_edge)
                edge_orientation[opposite_edge]=opposite_edge_orientation
                push!(δk,opposite_edge)
            else
                @assert edge_orientation[opposite_edge]==opposite_edge_orientation
            end
          end
      end
      δk_minus_one=copy(δk)
    end

end

function rotate_cell!(cell_vertices,
                      cell_vertices_coordinates,
                      cell_edges,
                      edge_vertices,
                      cell,
                      edge_orientation)
    starting_vertex_of_edge=Vector{Int}(undef,4)
    for l=1:4
      edge_gid=cell_edges[cell][l]
      @assert haskey(edge_orientation,edge_gid)
      if edge_orientation[edge_gid]==FORWARD
            starting_vertex_of_edge[l]=edge_vertices[edge_gid][1]
      else
            starting_vertex_of_edge[l]=edge_vertices[edge_gid][2]
      end
    end

    origin_vertex_of_cell=-1
    if ((starting_vertex_of_edge[3] == starting_vertex_of_edge[1]) ||
      (starting_vertex_of_edge[3] == starting_vertex_of_edge[2]))
      origin_vertex_of_cell = starting_vertex_of_edge[3]
    elseif ((starting_vertex_of_edge[4] == starting_vertex_of_edge[1]) ||
            (starting_vertex_of_edge[4] == starting_vertex_of_edge[2]))
      origin_vertex_of_cell = starting_vertex_of_edge[4]
    else
      @assert false
    end

    function rotate!(cell_wise_array,cell)
      tmp = cell_wise_array[cell][1];
      cell_wise_array[cell][1] = cell_wise_array[cell][2]
      cell_wise_array[cell][2] = cell_wise_array[cell][4]
      cell_wise_array[cell][4] = cell_wise_array[cell][3]
      cell_wise_array[cell][3] = tmp
    end

    while (cell_vertices[cell][1] != origin_vertex_of_cell)
      rotate!(cell_vertices,cell)
      rotate!(cell_vertices_coordinates,cell)
    end
end

edge_orientation=Dict{Int,Int8}()

next_cell_with_unoriented_edge=1

const starting_edges=[1,3]
const parallel_edges=[2,1,4,3]
const FORWARD=Int8(0)
const BACKWARD=Int8(1)

while (next_cell_with_unoriented_edge !=-1)
    for l in (1,2)
      lid_starting_edge=starting_edges[l]
      gid_starting_edge=cell_edges[next_cell_with_unoriented_edge][lid_starting_edge]
      if (!haskey(edge_orientation,gid_starting_edge))
          orient_one_set_of_parallel_edges!(edge_orientation,
                                          cell_vertices,
                                          cell_edges,
                                          edge_cells,
                                          edge_vertices_sorted,
                                          next_cell_with_unoriented_edge,
                                          lid_starting_edge)
      end
    end

    for l=1:length(cell_edges[next_cell_with_unoriented_edge])
      edge=cell_edges[next_cell_with_unoriented_edge][l]
      @assert haskey(edge_orientation,edge)
    end

    next_cell_with_unoriented_edge =
      get_next_unoriented_cell(edge_orientation,
                              cell_edges,
                              next_cell_with_unoriented_edge)
end

cell_vertices_new = deepcopy(collect(cell_vertices))
for cell=1:length(cell_vertices_new)
    rotate_cell!(cell_vertices_new,cell_vertices_coordinates,cell_edges,edge_vertices_sorted,cell,edge_orientation)
end



_topo = UnstructuredGridTopology(nodes_2d,Gridap.Arrays.Table(cell_vertices_new),cell_type,polytopes,Gridap.Geometry.NonOriented())

get_cell_permutations(_topo,1)

cell_vertices_new
