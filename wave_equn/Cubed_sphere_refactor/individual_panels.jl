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

function map_xy_2_XYZ(xy,panel)
  x,y = xy
  if panel == 1
    XYZ = Point(1.0,x,y)
  elseif panel == 2
    XYZ = Point(-x,1.0,y)
  elseif panel == 3
    XYZ = Point(-1.0,-x,y)
  elseif panel == 4
    XYZ = Point(x,-1.0,y)
  elseif panel == 5
    XYZ = Point(-y,x,1.0)
  elseif panel == 6
    XYZ = Point(y,x,-1.0)
  end
  map_cube_to_sphere(XYZ)
end

function map_cube_to_sphere(XYZ)
  x,y,z = XYZ
  xₛ = x*sqrt(1.0-y^2/2-z^2/2+y^2*z^2/(3.0))
  yₛ = y*sqrt(1.0-z^2/2-x^2/2+x^2*z^2/(3.0))
  zₛ = z*sqrt(1.0-x^2/2-y^2/2+x^2*y^2/(3.0))
  Point(xₛ,yₛ,zₛ)
end



function get_panel(panel,node_coordinates,cell_node_ids)

  mapped_node_coordinates = map( (x)-> map_xy_2_XYZ(x,panel), node_coordinates  )
  tab_cell_node_ids = Table(cell_node_ids)

  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_types = fill(1,length(tab_cell_node_ids))
  cell_reffes=[reffes]

  Gridap.Geometry.UnstructuredGrid(vec(mapped_node_coordinates),
                                  tab_cell_node_ids,
                                  cell_reffes,
                                  cell_types,
                                  Gridap.Geometry.NonOriented())

end



ncells_per_panel = 4
npanels = 6
ref_panel = CartesianDiscreteModel((-1,1,-1,1),(ncells_per_panel,ncells_per_panel))

cell_coordinates = get_cell_coordinates(ref_panel)
cell_node_ids = get_cell_node_ids(ref_panel)
node_coordinates = get_node_coordinates(ref_panel)

for panel = collect(1:npanels)
  cube_panel = get_panel(panel,node_coordinates,cell_node_ids)
  writevtk(cube_panel,datadir("CubedSphere")*"/panel$(panel)",append=false)
end
