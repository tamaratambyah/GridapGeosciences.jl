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
include("cube_surface_1_cell_per_panel.jl")
include("structs.jl")
include("refinement.jl")
include("maps.jl")

dir = datadir("CubedSphereRefactor/Test_problem")
!isdir(dir) && mkdir(dir)


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()
cube_model = UnstructuredDiscreteModel(cube_grid,topo,face_labels)
ref_cube_model = Gridap.Adaptivity.refine(cube_model)
cube_nodes = get_node_coordinates(cube_grid)

### Analytical map
CSgrid = CubedSphereGrid(cube_grid,map_cube_to_sphere)
_CSmodel = CubedSphereDiscreteModel(CSgrid,topo,face_labels)
CSmodel = Gridap.Adaptivity.refine(_CSmodel)

order = 1
Ω = Triangulation(CSmodel)
dΩ = Measure(Ω,2*order+1)
H1 = FESpace(CSmodel,ReferenceFE(lagrangian,VectorValue{3,Float64},order),
            conformity=:H1)

u0(x) = VectorValue(1,2,3)

u_cf = CellField(u0,Ω)
u_int = interpolate(u0,H1)

a(u,v) = ∫( u⋅v )dΩ
l(v) = ∫( v⋅u0  )dΩ
op = AffineFEOperator(a,l,H1,H1)
u_proj = solve(LUSolver(),op)

cf = ["u_cf"=>u_cf,"u_int"=>u_int,"u_proj"=>u_proj]
writevtk(Ω,dir*"/test",cellfields=cf,append=false)
