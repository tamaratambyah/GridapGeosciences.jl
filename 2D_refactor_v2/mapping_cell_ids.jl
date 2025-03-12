using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using BenchmarkTools
using LinearAlgebra
using FillArrays

include("initialise.jl")



# pick a model, and generate coordinates of panel 1
_model = ref_ref_ref_ref_model
panel_ids = get_panel_ids(_model)
cell_coords = get_cell_coordinates(_model)
p1 = findall(x->x==1,panel_ids)
cell_coords_panel1 = cell_coords[p1]
cell_coords_panel1_2D = lazy_map(BumpMap(), cell_coords_panel1)

num_cells_p1 = Int64(num_cells(_model)/6) #length(p1)

## create cartesian grid for panel 1
nC_p1 = Tuple(fill(Int(sqrt(num_cells_p1)),2))
panel1 = UnstructuredGrid(CartesianGrid((-1,1,-1,1), nC_p1 ))
_cell_coords_panel1 = get_cell_coordinates(panel1)

reindex = reorder_cell_ids(_model)
println(cell_coords_panel1_2D .== _cell_coords_panel1[reindex])




################################################################################
# reordering via recusion --- debugging
#############################################################################

# Dc = 2
# n_cells_per_panel = Int(num_cells(_model)/6)
# reindex = collect(1:n_cells_per_panel)



## now do the reordering for 1:4 partition rule
# n = Int(num_cells(_model.parent)/6)
# nc_fine = Tuple(fill(Int(sqrt(n)),Dc))
# f2c_cell_map, fcell_to_child_id = Gridap.Adaptivity._create_cartesian_f2c_maps(nc_fine, (2,2))
# c2f_cell_map = Adaptivity.get_o2n_faces_map(f2c_cell_map)
# p = collect(1: length(c2f_cell_map))

# reindex .= c2f_cell_map[p].data
# println(cell_coords_panel1_2D .== _cell_coords_panel1[reindex])

# og_cell_map = copy(c2f_cell_map)

# recurse
# prev_cell_map = copy(c2f_cell_map)


# n = Int(num_cells(_model.parent.parent)/6)
# nc_fine = Tuple(fill(Int(sqrt(n)),Dc))
# _f2c_cell_map,  = Gridap.Adaptivity._create_cartesian_f2c_maps(nc_fine, (2,2))
# _c2f_cell_map = Adaptivity.get_o2n_faces_map(_f2c_cell_map)

# p = sortperm(_f2c_cell_map)
# reindex .= c2f_cell_map[p].data
# println(cell_coords_panel1_2D .== _cell_coords_panel1[reindex])


# p = sortperm(_f2c_cell_map)
# o_reindex = prev_cell_map[p].data
# if length(p) == length(og_cell_map)
#   reindex .= og_cell_map[p].data
# else
#   reindex .= og_cell_map[o_reindex].data
# end

# println(cell_coords_panel1_2D .== _cell_coords_panel1[reindex])


# recurse again
# prev_cell_map = copy(_c2f_cell_map)


# n = Int(num_cells(_model.parent.parent.parent)/6)
# nc_fine = Tuple(fill(Int(sqrt(n)),Dc))
# o_f2c_cell_map,  = Gridap.Adaptivity._create_cartesian_f2c_maps(nc_fine, (2,2))
# o_c2f_cell_map = Adaptivity.get_o2n_faces_map(o_f2c_cell_map)

# p = sortperm(o_f2c_cell_map)
# o_reindex = prev_cell_map[p].data

# if length(p) == length(og_cell_map)
#   reindex .= og_cell_map[p].data
# else
#   reindex .= og_cell_map[o_reindex].data
# end


# println(cell_coords_panel1_2D .== _cell_coords_panel1[reindex])
