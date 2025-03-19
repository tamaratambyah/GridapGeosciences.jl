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
include("maps/cube_cell_map.jl")
include("maps/sphere_cell_map.jl")
include("adaptivity_tools/panel_ids_from_refinement.jl")
include("cube_topo/cube_surface_1_cell_per_panel.jl")
include("helpers.jl")
include("ManifoldGrid.jl")

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

a = π/4.0
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
rp1 = map(TensorValue,rotate_panel_p_to_1)
r1p = map(TensorValue,rotate_panel_1_to_p)

_A , _B, _b = bump_matrics(a)
A_bump = TensorValue(_A)
B_bump = TensorValue(_B)
b_bump = VectorValue(_b)

PanelMap() = PanelRotationMap(rp1)
InvPanelMap() = PanelRotationMap(r1p)
BumpMap() = Panel1BumpMap(A_bump,B_bump,b_bump)
Sigma() = SigmaMap(r)

Rp1 = map(p->PanelRotationMap(rp1[p]), 1:6)
R1p = map(p->PanelRotationMap(r1p[p]), 1:6)
Bump = Panel1BumpMap(A_bump,B_bump,b_bump)
Gnomonic = GnomonicMap()
SSigma = SigmaMap(r)


CubePhysCellMap() = CubePhysCellMap(Rp1,R1p,Bump)
SphereAmbientCellMap() = SphereAmbientCellMap(Rp1,R1p,Bump,Gnomonic,SSigma)
