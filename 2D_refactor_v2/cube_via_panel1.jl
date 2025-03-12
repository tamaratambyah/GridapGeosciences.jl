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
function cube_nodes(model)
  panel_ids = get_panel_ids(model)

  reindex = reorder_cell_ids(model)
  num_cells_p1 = Int64(num_cells(model)/6)

  ## create cartesian grid for panel 1 and reindex
  nC_p1 = Tuple(fill(Int(sqrt(num_cells_p1)),2))

  domain = (-1.0,1.0,-1.0,1.0) # atan(-1), atan(1)
  panel1 = UnstructuredGrid(CartesianGrid(domain, nC_p1 ))
  _local_coords = get_cell_coordinates(panel1)

  local_coords = _local_coords[reindex]

  ids = repeat(Base.OneTo(num_cells_p1),6)
  all_local_coords = lazy_map(Reindex(local_coords),ids)


  # return coords
  _coords_3D = lazy_map(BumpMap(), all_local_coords)
  coords_3D = lazy_map(InvPanelMap(), _coords_3D, panel_ids)

  return all_local_coords,  coords_3D

end


# pick a model, and generate coordinates of panel 1
_model = ref_ref_model
all_local_coords,  coords_3D, = cube_nodes(_model)

test_coords(coords_3D,get_cell_coordinates(_model))
