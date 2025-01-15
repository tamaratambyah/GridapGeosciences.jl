using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Test
using LinearAlgebra
using FillArrays
# include("individual_panels.jl")

function map_xy_2_XYZ_ccam(xy,panel)
  x,y = xy
  if panel == 0
    XYZ = Point(1.0, x, y)
  elseif panel == 1
    XYZ = Point(-y, x, 1.0)
  elseif panel == 2
    XYZ = Point(-y, 1.0, -x)
  elseif panel == 3
    XYZ = Point(-1.0, -y, -x)
  elseif panel == 4
    XYZ = Point(x, -y, -1.0)
  elseif panel == 5
    XYZ = Point(x, -1.0, y)
  end
  XYZ
  # map_cube_to_sphere(XYZ)
end

function get_face_ids(corner_idx)
  south = collect(corner_idx[1]:corner_idx[2])
  north = collect(corner_idx[3]:corner_idx[4])
  west = collect(corner_idx[1]:3:corner_idx[3])
  east = collect(corner_idx[2]:3:corner_idx[4])
  return south,west,north,east
end


n = 2

ref_panel = UnstructuredGrid( get_grid( CartesianDiscreteModel((-1,1,-1,1),(n,n)) ))
ref_grid = get_grid( CartesianDiscreteModel((-1,1,-1,1),(n,n)))

# panel 0
xa = get_node_coordinates(ref_panel)
Nx = Int32(length(xa))
nx = Int32(sqrt(Nx))
corner_idx0 = [1, n, Nx-nx+1, Nx] # bottom left, bottom right, top left, top right

# join panels 0 and 1
#################################################################################
p = 0
ids0 = get_cell_node_ids(ref_panel)
nodes0 = get_node_coordinates(ref_panel)
mapped_nodes0 = map( (xy)-> map_xy_2_XYZ_ccam(xy,p), nodes0  )

p = 1
ids1 = Table( map(Broadcasting(i -> (i+Nx) - nx), ids0) )
mapped_nodes1 = map( (xy)-> map_xy_2_XYZ_ccam(xy,p), nodes0  )
nodes_to_append = mapped_nodes1[nx+1:Nx]


new_nodes = lazy_append(mapped_nodes0,nodes_to_append)./1
new_cell_ids = append_tables_globally(ids0,ids1)


reffes = LagrangianRefFE(Float64,QUAD,1)
cell_types = fill(1,length(new_cell_ids))
cell_reffes=[reffes]

append_grid = Gridap.Geometry.UnstructuredGrid(new_nodes,
                                new_cell_ids,
                                cell_reffes,
                                cell_types,
                                Gridap.Geometry.NonOriented())

writevtk(append_grid,datadir("CubedSphere")*"/append_grid",append=false)

#################################################################################





#################################################################################
p = 2
top_right2 = maximum(ids1.data) + Nx -(2*n-1)
dif = top_right2-top_right0
ids2 = Table( map(Broadcasting(i -> (i+dif)), ids0) )


p = 3
ids3 = map(Broadcasting(i -> (i+Nx)-2*n + 1), ids2)

p = 4
ids4 = map(Broadcasting(i -> (i+Nx)-(2*n-1)-1), ids3)


# south_ids,west_ids,north_ids,east_ids = get_face_ids(corner_idx0)


# data = idsa.data
# _data = map( x -> (x == 7 ? north_ids1[1] : x), data)

# for (i,cell_list) in enumerate(idsa)
#   println(cell_list)
#   for k in 1:length(cell_list)
#     coord = cell_list[k]

#     for idx in 1:length(north_ids)
#       if coord == north_ids[idx]
#         println("replacing")
#         idsa[i][k] = south_ids1[idx ]
#       end
#     end

#     for idx in 1:length(east_ids)
#       if coord == east_ids[idx]
#         idsa[i][k] = west_ids1[idx ]
#       end
#     end

#   end
# end
# corner_idx1 = map(idx->idx + p*Nx, corner_idx0)
# south_ids1,west_ids1,north_ids1,east_ids1 = get_face_ids(corner_idx1)
