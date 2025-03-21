using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Fields
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

include("initialise.jl")

_model = ref_model
# _model = ref_model
manifold_grid = CubeGrid(_model)
panel_ids = get_panel_ids(manifold_grid)

order = 4

Ω = BodyFittedTriangulation(_model, manifold_grid, IdentityVector(num_cells(manifold_grid)))
# Ω = GenericTriangulation(manifold_grid)
dΩ = Measure(Ω,order)

f = CellField(1,Ω)

quad = CellQuadrature(Ω,order) #dΩ.quad
trian_f = get_triangulation(f)
trian_x = get_triangulation(quad)

@Gridap.Helpers.check is_change_possible(trian_f,trian_x)

b = change_domain(f,quad.trian,quad.data_domain_style)
x = get_cell_points(quad)
bx = b(x)


# cell_map = get_cell_map(manifold_grid)
# cell_Jt = lazy_map(∇,cell_map)
# cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)


grid = manifold_grid.topo_grid
cell_to_shapefuns = get_cell_shapefuns(grid)

cell_to_coords = get_cell_coordinates(grid)
_mapp_coords = lazy_map(PanelMap(),cell_to_coords,panel_ids)
mapp_coords = lazy_map(BumpMap(),_mapp_coords)



latlon_panel1 = lazy_map(GnomonicMap(), mapp_coords)
sphere_panel1 = lazy_map(Sigma(),latlon_panel1)
sphere_panelp = lazy_map(InvPanelMap(), sphere_panel1, panel_ids)


_cell_map = lazy_map(linear_combination,sphere_panelp,cell_to_shapefuns)

_cell_Jt = lazy_map(∇,_cell_map)
_cell_Jtx = lazy_map(evaluate,_cell_Jt,quad.cell_point)
int = lazy_map(IntegrationMap(),bx,quad.cell_weight,_cell_Jtx)

sum(int)
# 6*(2*a)^2
4*π*r^2
