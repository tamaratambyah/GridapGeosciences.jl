"""
Basic test of the generic MappedDiscreteModel in Gridap.
Take a 2D point and map it to another 2D point via stretch or compression
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

map(x) = 2*x #sin(x[1]...)  #x⋅(1-x)

p = 1
u_ex(x) = VectorValue(2*x[1],2*x[2])

domain = (0,4,0,4)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
uh = fe_problem_Hdiv(model,p,u_ex)
writevtk(Triangulation(model),datadir("Mapped")*"/2D_linear_og",cellfields=["u"=>uh],append=false)

new_model = MappedDiscreteModel(model,map)
uh_map = fe_problem_Hdiv(new_model,p,u_ex)
writevtk(Triangulation(new_model),datadir("Mapped")*"/2D_linear_mapped",cellfields=["u_map"=>uh_map],append=false)


## compare values
pt = Point(1,1)

grid = get_grid(model)
nodes = println(get_node_coordinates(grid))
dofs = println(get_free_dof_values(uh))
cmaps = get_cell_map(model)#./1
evaluate(cmaps[4],pt)

cmaps = Gridap.Arrays.collect1d(get_cell_map(Triangulation(model)))
cell_Jt = Gridap.Arrays.collect1d( lazy_map(∇,cmaps))
cell_Jtx = evaluate(cell_Jt[1],pt)


new_grid = get_grid(new_model) #MappedGrid(grid,map)
new_nodes = println(get_node_coordinates(new_grid))
new_dofs = println(get_free_dof_values(uh_map))
new_cmaps = get_cell_map(new_model)./1
evaluate(new_cmaps[2],pt)

new_cmaps = Gridap.Arrays.collect1d(get_cell_map(Triangulation(new_model)))
new_cell_Jt = Gridap.Arrays.collect1d( lazy_map(∇,new_cmaps))
new_cell_Jtx = evaluate(new_cell_Jt[2],pt)
