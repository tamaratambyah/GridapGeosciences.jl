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
include("maps/map_matrices.jl")
include("maps/maps.jl")
include("adaptivity_tools/panel_ids_from_refinement.jl")
include("cube_topo/cube_surface_1_cell_per_panel.jl")
# include("helpers.jl")

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

a = π/4
r = a*sqrt(3)

cube_model_3D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel(a)...)

##### make some refined models
model = cube_model_3D
ref_model = Gridap.Adaptivity.refine(model)
ref_ref_model = Gridap.Adaptivity.refine(ref_model)
ref_ref_ref_model = Gridap.Adaptivity.refine(ref_ref_model)
ref_ref_ref_ref_model = Gridap.Adaptivity.refine(ref_ref_ref_model)
ref_ref_ref_ref_ref_model = Gridap.Adaptivity.refine(ref_ref_ref_ref_model)


## get rotation maps
rotate_panel_p_to_1, rotate_panel_1_to_p = panel_rotations()
Rp1 = map(TensorValue,rotate_panel_p_to_1)
R1p = map(TensorValue,rotate_panel_1_to_p)

_A , _B, _b = bump_matrics(a)
A_bump = TensorValue(_A)
B_bump = TensorValue(_B)
b_bump = VectorValue(_b)

PanelMap() = PanelRotationMap(Rp1)
InvPanelMap() = PanelRotationMap(R1p)
BumpMap() = Panel1BumpMap(A_bump,B_bump,b_bump)


Sigma() = SigmaMap(r)
