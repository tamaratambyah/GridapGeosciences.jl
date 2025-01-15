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

# function is_upper(coords)

#   for k in coords
#     x,y = k
#     if x == -1 || y == -1
#       return false
#     else
#       return true
#     end
#   end
# end




n = 2

ref_panel =  CartesianDiscreteModel((-1,1,-1,1),(n,n))
ref_grid =  get_grid( ref_panel )
ref_topo = get_grid_topology( ref_panel  )


labels = get_face_labeling(ref_panel)
boundary_ids = get_face_mask(labels,"boundary",0)
mask = findall(boundary_ids)

faces = Gridap.Geometry.get_faces(ref_topo,2,2)

cell_node_ids = get_cell_node_ids(ref_panel)
cell_coordinates = get_cell_coordinates(ref_panel)

mask = map(east,cell_coordinates)

function east(coords)

  for k in coords
    x,y = k
    if x == 1
      return true
    else
      return false
    end
  end
end


edge_grid = GridPortion(ref_grid,mask)
