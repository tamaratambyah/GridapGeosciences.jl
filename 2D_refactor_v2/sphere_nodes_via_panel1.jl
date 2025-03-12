using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools
using Plots

include("initialise.jl")


## above as a function:
function cubed_sphere_nodes(model)
  panel_ids = get_panel_ids(model)

  reindex = reorder_cell_ids(model)
  num_cells_p1 = Int64(num_cells(model)/6)

  ## create cartesian grid for panel 1 and reindex
  nC_p1 = Tuple(fill(Int(sqrt(num_cells_p1)),2))

  cangels_domain = (-π/4,π/4,-π/4,π/4) # atan(-1), atan(1)
  panel1 = UnstructuredGrid(CartesianGrid(cangels_domain, nC_p1 ))
  _central_angels_coords = get_cell_coordinates(panel1)

  central_angels_coords = _central_angels_coords[reindex]

  ids = repeat(Base.OneTo(num_cells_p1),6)
  all_central_angels_coords = lazy_map(Reindex(central_angels_coords),ids)


  # return local cartesian coords (non-uniform)
  local_coords = lazy_map(InverseCentralAngleMap(), all_central_angels_coords)


  # return lat lons
  panel1_latlon = lazy_map(OtherGnomonicMap(), all_central_angels_coords)
  panel1_sphere = lazy_map(Sigma(),panel1_latlon)
  panelp_sphere = lazy_map(InvPanelMap(), panel1_sphere, panel_ids)
  panelp_latlon = lazy_map(Sigma(), panelp_sphere)


  return all_central_angels_coords,  panelp_sphere, local_coords, panelp_latlon

end


# pick a model, and generate coordinates of panel 1
_model = model
all_central_angels_coords,  panelp_sphere, local_coords, panelp_latlon = cubed_sphere_nodes(_model)


###############################################################################
### DEBUGGING


# panel_ids = get_panel_ids(_model)
# cell_coords = get_cell_coordinates(_model)
# p1 = findall(x->x==1,panel_ids)
# cell_coords_panel1 = cell_coords[p1]
# cell_coords_panel1_2D = lazy_map(BumpMap(), cell_coords_panel1)

# num_cells_p1 = Int64(num_cells(_model)/6) #length(p1)

# ## create cartesian grid for panel 1 and reindex
# nC_p1 = Tuple(fill(Int(sqrt(num_cells_p1)),2))
# panel1 = UnstructuredGrid(CartesianGrid((-1,1,-1,1), nC_p1 ))
# _cell_coords_panel1 = get_cell_coordinates(panel1)

# reindex = reorder_cell_ids(_model)
# println(cell_coords_panel1_2D .== _cell_coords_panel1[reindex])

# ids = repeat(Base.OneTo(num_cells_p1),6)
# all_panel1 = lazy_map(Reindex(_cell_coords_panel1[reindex]),ids)

# coords_panel1_3D = lazy_map(BumpMap(), all_panel1)
# coords_panelp = lazy_map(InvPanelMap(), coords_panel1_3D, panel_ids)
# test_coords(coords_panelp,cell_coords)


# _reindex() = reorder_cell_ids(_model)
# @benchmark _reindex()

# cache = array_cache(all_panel1)
# bm() = lazy_collect(cache,all_panel1)
# @benchmark bm()



# # now do the above for central angles
# cangels_domain = (-π/4,π/4,-π/4,π/4) # atan(-1), atan(1)
# panel1 = UnstructuredGrid(CartesianGrid(cangels_domain, nC_p1 ))
# _central_angels_coords = get_cell_coordinates(panel1)
# central_angels_coords = _central_angels_coords[reindex]

# panel1_ids = ones(Int,size(p1))
# plot_coords(central_angels_coords,panel1_ids,"panel1_central_angles","alpha","beta")

# # return local cartesian coords (non-uniform)
# _coords_panel1_2D = lazy_map(InverseCentralAngleMap(), central_angels_coords)
# cache = array_cache(_coords_panel1_2D)
# bm1() = lazy_collect(cache,_coords_panel1_2D)
# @benchmark bm1()

# plot_coords(_coords_panel1_2D,panel1_ids,"panel1_local_coords","x","y")

# # return lat lons
# panel1_latlon = lazy_map(OtherGnomonicMap(), central_angels_coords)

# lazy_map(x->panel1_latlon,1:6)./1

# B = lazy_append(panel1_latlon,panel1_latlon)
# C = lazy_append(B,B)
# D = lazy_append(C,B)
# cache = array_cache(D)
# bm() = lazy_collect(cache,D)
# @benchmark bm()



# all_panel1_latlon = lazy_append(C,B) #lazy_map(x->panel1_latlon,1:6) # repeat(panel1_latlon,6)
# panel1_sphere = lazy_map(Sigma(),all_panel1_latlon)
# panelp_sphere = lazy_map(InvPanelMap(), panel1_sphere, panel_ids)
# panelp_latlon = lazy_map(Sigma(), panelp_sphere)
# cache = array_cache(panelp_latlon)
# bm2() = lazy_collect(cache,panelp_latlon)
# @benchmark bm2()

# plot_coords(panel1_latlon,panel1_ids,"panel1_latlon","longitude","latitude") #### plot the lat lons in panel 1
# plot_coords(panelp_latlon,panel_ids,"latlon","longitude","latitude") #### plot the lat lons in panel i=1,..,6
# plot_coords(panelp_sphere,panel_ids,"sphere","longitude","latitude") #### plot the lat lons in panel i=1,..,6


# mapped_grid = make_grid(get_grid_topology(_model),panelp_sphere)
# writevtk(mapped_grid,dir*"/cubed_sphere",append=false)
