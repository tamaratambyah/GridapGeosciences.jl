using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../../../src/initialise.jl")

P1 = one(TensorValue{2,2,Float64})
P2 = TensorValue{2,2,Float64}(-1.0, 0.0, 0.0, 1.0)
Ps = [P1;P2]
Rs = [rp1_3D[1];rp1_3D[5]]

_P2 = TensorValue{2,2,Float64}(-1.0, 0.0, 0.0, -1.0)
_Ps = [P1,_P2]

shifts = [VectorValue(0.0,0.0);VectorValue(0,π/2)]
invshifts = [VectorValue(0.0,0.0);VectorValue(0,-π/2)]

######### test panels 1 & 2
RADIUS = sqrt(3)
# _model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π/4,π/4, -3*π/4,π/4),(8,16)))
_model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π/4,π/4, -3*π/4,π/4),(4,8)))
# _model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π/4,π/4, -3*π/4,π/4),(2,4)))
# _model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π/4,π/4, -3*π/4,π/4),(1,2)))

# panel_ids = [fill(2,64);fill(1,64)]
panel_ids = [fill(2,16);fill(1,16)]
# panel_ids = [fill(2,4);fill(1,4)]
# panel_ids = [fill(2,1);fill(1,1)]

grid = get_grid(_model)
cmaps = get_cell_map(grid)


k = lazy_map(p-> InversionField(Ps[p]) ∘ ShiftField(shifts[p]) ∘ InversionField(_Ps[p]), panel_ids)
_cmaps = lazy_map(∘,k,cmaps)

_grid = UnstructuredGrid(get_node_coordinates(grid),get_cell_node_ids(grid),get_reffes(grid),get_cell_type(grid),
  OrientationStyle(grid),nothing,_cmaps)

model = UnstructuredDiscreteModel(_grid,get_grid_topology(_model),get_face_labeling(_model))

######### FE
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model;tags="boundary")
nΓ = get_normal_vector(Γ)
writevtk(Ω,dir*"/test",append=false)

dΩ = Measure(Ω,10)
dΓ = Measure(Γ,10)

V = TestFESpace(model, ReferenceFE(lagrangian,Float64,3); conformity=:H1)
U = TrialFESpace(V)

uX_scalar(x) = x[1]*x[2]*x[3]

function u_scalar_ambient2parametric2(p::Int,uX::Function)
  function _u(αβ)
    θϕ = GnomonicField()(αβ)
    _XYZ = SigmaField(RADIUS)(θϕ)
    XYZ = PanelRotationField(Rs[p])(_XYZ)
    uX(XYZ)
  end
end

cell_field = map(p->GenericField(u_scalar_ambient2parametric2(p,uX_scalar)),panel_ids)
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())
writevtk(Ω,dir*"/test_u",cellfields=["u"=>ucf],append=false)


M = CellField(metric1,Ω)
invM = CellField(invmetric1,Ω)
sqrtM = CellField(sqrtmet1,Ω)

mNew =  Metric(M,sqrtM,invM,metric1,sqrtmet1,invmetric1)

lapucf = surface_laplacian(ucf,mNew) #1/sqrtM * divergence( sqrtM*( invM ⋅ gradient(ucf) ) )

rhs = ucf + lapucf
sum(∫(rhs)dΩ)

mass(u,v) = ∫( (u*v)*sqrtM )dΩ
stiffnes(u,v) = ∫( (∇(v)⋅ (invM⋅ ∇(u) )) *sqrtM)dΩ
bound(v) = ∫( (((invM ⋅∇(ucf))⋅(nΓ)*v)*sqrtM ) )dΓ
force(v) = ∫(  (rhs*v)*sqrtM )dΩ

helmholtz_biform(u,v) = mass(u,v) - stiffnes(u,v)
helmholtz_liform(v) = force(v) - bound(v)

op = AffineFEOperator(helmholtz_biform,helmholtz_liform,U,V)

uh = solve(LUSolver(),op)

e = l2(uh-ucf,dΩ)

mask = (panel_ids.==2)
Ωp = Triangulation(model,mask)
writevtk(Ωp,dir*"/test_u",cellfields=["u"=>ucf,"uh"=>uh,"e"=>ucf-uh,],append=false)

######## map to proper parametric space

invMapping = map(p->InversionField(_Ps[p]) ∘  ShiftField(invshifts[p]) ∘ InversionField(Ps[p]) , panel_ids)

_Ω = Triangulation(_model)
_pts = get_cell_points(_Ω)

_cf  = change_domain(uh.cell_field,ReferenceDomain(),PhysicalDomain())
cf_mapped = lazy_map(Broadcasting(∘),get_data(_cf),invMapping)
uh_mapped = CellData.GenericCellField(cf_mapped,_Ω,PhysicalDomain() )
uh_mapped(_pts)

_cf_mapped = lazy_map(Broadcasting(∘),get_data(ucf),invMapping)
ucf_mapped = CellData.GenericCellField(_cf_mapped,_Ω,PhysicalDomain() )

ucf_mapped(_pts)

function _u_scalar_ambient2parametric2(p::Int,uX::Function)
  function _u(αβ)
    αβ1 = ShiftField(shifts[p])(αβ )
    θϕ = GnomonicField()(αβ1)
    _XYZ = SigmaField(RADIUS)(θϕ)
    XYZ = PanelRotationField(r1p_3D[p])(_XYZ)
    uX(XYZ)
  end
end

_rot_cf = map(p->GenericField(_u_scalar_ambient2parametric2(p,uX_scalar)),panel_ids)
rot_ucf = CellData.GenericCellField(_rot_cf,_Ω,PhysicalDomain())

writevtk(Triangulation(_model),dir*"/test_u_mapped",
cellfields=["u"=>ucf_mapped,"uh"=>uh_mapped,"e"=>ucf_mapped-uh_mapped,"urot"=>rot_ucf,"eu"=>ucf_mapped-rot_ucf],append=false)
