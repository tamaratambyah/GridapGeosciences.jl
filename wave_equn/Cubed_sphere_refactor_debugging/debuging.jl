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
  x,y=xy
  if panel==1
    XYZ=Point(1.0,x,y)
  elseif panel==2
    XYZ=Point(-x,1.0,y)
  elseif panel==3
    XYZ=Point(-1.0,-x,y)
  elseif panel==4
    XYZ=Point(x,-1.0,y)
  elseif panel==5
    XYZ=Point(-y,x,1.0)
  elseif panel==6
    XYZ=Point(y,x,-1.0)
  end
  XYZ
end

function generate_ptr(n)
  nvertices = 4
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end

ncells_per_panel = 4
npanels = 6
pol = QUAD

ref_panel = CartesianDiscreteModel((-1,1,-1,1),(ncells_per_panel,ncells_per_panel))

cell_coordinates = get_cell_coordinates(ref_panel)
cell_node_ids = get_cell_node_ids(ref_panel)
node_coordinates = get_node_coordinates(ref_panel)

panel = 1

# mapped_node_coordinates = map( (x)-> x, node_coordinates  )
  mapped_node_coordinates = map( (x)-> map_xy_2_XYZ(x,panel), node_coordinates  )
  cell_vertex_lids = Table(cell_node_ids)


  reffes = LagrangianRefFE(Float64,pol,1)
  cell_types = fill(1,length(cell_vertex_lids))
  cell_reffes=[reffes]
  cube_face_grid = Gridap.Geometry.UnstructuredGrid(vec(mapped_node_coordinates),
                                          cell_vertex_lids,
                                          cell_reffes,
                                          cell_types,
                                          Gridap.Geometry.NonOriented())



  writevtk(cube_face_grid,datadir("CubedSphere")*"/panel$(panel)",append=false)

# endd



########### panels
data = [ 2,4,6,8, 4,3,8,7, 3,1,7,5, 1,2,5,6, 6,8,5,7, 1,3,2,4  ]
cube_vertex_ids = Gridap.Arrays.Table(data,ptr)


node_coordinates = Vector{Point{2,Float64}}(undef,8)
for i in 1:length(node_coordinates)
  node_coordinates[i]=Point{2,Float64}(0.0,0.0)
end

reffes = LagrangianRefFE(Float64,pol,1)
cell_types = fill(1,length(cube_vertex_ids))
cell_reffes=[reffes]
grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                        cube_vertex_ids,
                                        cell_reffes,
                                        cell_types,
                                        Gridap.Geometry.NonOriented())


model = Gridap.Geometry.UnstructuredDiscreteModel(grid)
writevtk(Triangulation(model),datadir("CubedSphere")*"/panels",append=false)



cube  = CartesianDiscreteModel((-1,1,-1,1,-1,1),(1,1,1))
writevtk(cube,datadir("CubedSphere")*"/cube",append=false)





##### OLD -- make nodes via cell coordinates

# lazy_map(Reindex(node_coordinates),mapped_node_coordinates)

# _mapped_node_coordinates = AbstractArray{Point{3,Float64},3}

# T = VectorValue{2,Float64}
# N = length(cell_coordinates)
# mapped_cell_coordinates  = fill(Vector{T}(undef,4),N)

# # map(cell_coordinates) do (cell_coords)
#   for i = 1:N
#     mapped_cell_coordinates[i] = map( (x)-> x, cell_coordinates[i])
#     # mapped_cell_coordinates[i] = map( (x)-> map_xy_2_XYZ(x,panel), cell_coords)
#   end
# # end

# new_node_coordinates = Vector{T}(undef,length(node_coordinates))
# Gridap.Geometry._cell_vector_to_dof_vector!(new_node_coordinates,cell_node_ids,mapped_cell_coordinates)


# mapped_node_coordinates = Vector{T}(undef,length(node_coordinates))
# Gridap.Geometry._cell_vector_to_dof_vector!(mapped_node_coordinates,cell_node_ids,mapped_cell_coordinates)



# function generate_ptr(n)
#   nvertices = 4
#   ptr  = Vector{Int}(undef,n+1)
#   ptr[1]=1
#   for i=1:n
#     ptr[i+1]=ptr[i]+nvertices
#   end
#   ptr
# end

# ptr = generate_ptr(npanels)
# data = [ 2,4,6,8, 4,3,8,7, 3,1,7,5, 1,2,5,6, 6,8,5,7, 1,3,2,4  ]
# cube_vertex_ids = Gridap.Arrays.Table(data,ptr)


# topo = UnstructuredGridTopology(panel1)
