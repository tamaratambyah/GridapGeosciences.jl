"""
Take a 1D point and map it to a 2D point
The 1D point can be thought of as the parametrising coordinate.
i.e. a circle can be parameterised by a single coordinate, θ, as [sin(θ),cos(θ)]
"""

using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Test

include("feproblem.jl")
include("structs.jl")

# map(x) = VectorValue(x...,4*x...)
map(x) = VectorValue(cos(x[1]...),sin(x[1]...))
evaluate(map,Point(1.0))

p = 1
u_ex(x) = 2*x[1]

domain = (0,2*π,)
partition = (10)
# model = CartesianDiscreteModel(domain,partition)
model = CartesianDiscreteModel(domain,partition,isperiodic=(true,))
grid = get_grid(model)

uh = fe_problem(model,p,u_ex)
writevtk(Triangulation(model),datadir("Mapped")*"/1D_curve_og",cellfields=["u"=>uh],append=false)


mapp_grid = MyMappedGrid(grid,map)
new_model = MyMappedDiscreteModel(model,mapp_grid)
uh_map = fe_problem(new_model,p,u_ex)
writevtk(Triangulation(new_model),datadir("Mapped")*"/1D_curve_mapped",cellfields=["u_map"=>uh_map],append=false)


## compare values
pt = Point(0)

grid = get_grid(model)
nodes = get_node_coordinates(grid)
dofs = get_cell_dof_values(uh)
cmaps = get_cell_map(model)#./1
evaluate(cmaps[1],pt)

cmaps = Gridap.Arrays.collect1d(get_cell_map(Triangulation(model)))
cell_Jt = Gridap.Arrays.collect1d( lazy_map(∇,cmaps))
cell_Jtx = evaluate(cell_Jt[2],pt)

pt = Point(1)
new_grid = get_grid(new_model) #MappedGrid(grid,map)
new_nodes = get_node_coordinates(new_grid)
new_dofs = get_free_dof_values(uh_map)
new_cmaps = get_cell_map(new_model)./1
evaluate(new_cmaps[2],pt)

new_cmaps = Gridap.Arrays.collect1d(get_cell_map(Triangulation(new_model)))
new_cell_Jt = Gridap.Arrays.collect1d( lazy_map(∇,new_cmaps))
new_cell_Jtx = evaluate(new_cell_Jt[2],pt)







#### low level construction #####
cell_node_ids = get_cell_node_ids(grid)
old_nodes = get_node_coordinates(grid)

dims = num_dims(grid) + 1
T = VectorValue{dims,Float64} # eltype(old_nodes)

node_coordinates = Vector{T}(undef,length(old_nodes))
c_coor = get_cell_coordinates(grid)./1

N = length(c_coor)
map_c_coor  = fill(Vector{T}(undef,length(c_coor)),N)
for i = 1:N
  # map_c_coor[i] = Base.map( (x)->evaluate(map,x), c_coor[i])
  map_c_coor[i] = Base.map( (x)->map(x), c_coor[i])
end
println(map_c_coor)
println(c_coor)
g = lazy_map(  (x)-> map(x), c_coor)


Gridap.Geometry._cell_vector_to_dof_vector!(node_coordinates,cell_node_ids,map_c_coor)
node_coords = collect(node_coordinates)

model_map=get_cell_map(grid)
# geo_map=lazy_map(∘,map,model_map)

#### structs #####


# function Gridap.Visualization._vtkpoints(trian)
#   D = num_point_dims(trian)
#   if D == 1
#     D += 1
#   end
#   println(D)
#   x = get_node_coordinates(trian)
#   xflat = collect(x)
#   reshape(reinterpret(Float64,xflat),(D,length(x)))
# end

# writevtk(new_model,datadir("Mapped")*"/1D_curve",append=false)
