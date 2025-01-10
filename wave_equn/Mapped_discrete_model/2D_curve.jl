"""
Take a 2D point and map it to a 3D point
The 2D point can be thought of as the parametrising coordinate.
i.e. a sphere can be parameterised by [x,y]↦[ cos(x)cos(y), cos(x)sin(y), sin(x) ]
for x ∈ [-π/2, π/2 ] and y ∈ [0, 2π]
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

# map(x) = VectorValue(x[1],4*x[1],10*x[2])
map(x) = VectorValue(cos(x[1])*cos(x[2]),sin(x[1])*cos(x[2]),sin(x[2]))
evaluate(map,Point(1.0,2.0))

p = 1
u_ex(x) = 2*x[1]

domain = (0, 2*π,-π/2, π/2)
partition = (20,20)
model = CartesianDiscreteModel(domain,partition,isperiodic=(true,true))
writevtk(model,datadir("Mapped")*"/2D_grid",append=false)
grid = get_grid(model)
uh = fe_problem(model,p,u_ex)
writevtk(Triangulation(model),datadir("Mapped")*"/2D_curve_og",cellfields=["u"=>uh],append=false)


mapp_grid = MyMappedGrid(grid,map)
new_model = MyMappedDiscreteModel(model,mapp_grid)
writevtk(new_model,datadir("Mapped")*"/2D_curve",append=false)
uh_map = fe_problem(new_model,p,u_ex)
writevtk(Triangulation(new_model),datadir("Mapped")*"/2D_curve_mapped",cellfields=["u_map"=>uh_map],append=false)
