using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../../../src/initialise.jl")

######### test panels 1 & 2
RADIUS = sqrt(3)
_model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π/4,π/4, -π/4,3*π/4),(8,16)))
# _model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π/4,π/4, -π/4,3*π/4),(4,8)))
# _model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π/4,π/4, -π/4,3*π/4),(2,4)))
# _model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π/4,π/4, -π/4,3*π/4),(1,2)))

panel_ids = [fill(1,64);fill(2,64)]
# panel_ids = [fill(1,16);fill(2,16)]
# panel_ids = [fill(1,4);fill(2,4)]
# panel_ids = [fill(1,1);fill(2,1)]

grid = get_grid(_model)
cmaps = get_cell_map(grid)
evaluate(cmaps[6],Point(1,1))

k = lazy_map(p-> InversionField(Ps[p]) ∘ ShiftField(p), panel_ids)
_cmaps = lazy_map(∘,k,cmaps)
evaluate(_cmaps[8],Point(1,0))

_grid = UnstructuredGrid(get_node_coordinates(grid),get_cell_node_ids(grid),get_reffes(grid),get_cell_type(grid),
  OrientationStyle(grid),nothing,_cmaps)

model = UnstructuredDiscreteModel(_grid,get_grid_topology(_model),get_face_labeling(_model))


Ω = Triangulation(model)
Γ = BoundaryTriangulation(model;tags="boundary")
nΓ = get_normal_vector(Γ)
writevtk(Ω,dir*"/test",append=false)

dΩ = Measure(Ω,10)
dΓ = Measure(Γ,10)

V = TestFESpace(model, ReferenceFE(lagrangian,Float64,3); conformity=:H1)
U = TrialFESpace(V)

uX_scalar(x) = x[1]*x[2]*x[3]

P1 = one(TensorValue{2,2,Float64})
P2 = TensorValue{2,2,Float64}(1,0.0,0.0,-1)
Ps = [P1;P2]

function u_scalar_ambient2parametric2(p::Int,uX::Function)
  function _u(αβ)
    # αβ1 =  ShiftField(p)(αβ)
    # _αβ1 = InversionField(Ps[p])(αβ1)
    θϕ = GnomonicField()(αβ)
    # _θϕ = InversionField(Ps[p])(θϕ)
    _XYZ = SigmaField(RADIUS)(θϕ)
    # XYZ = PanelRotationField(r1p_3D[p])(_XYZ)

    uX(_XYZ)
  end
end

cell_field = map(p->GenericField(u_scalar_ambient2parametric2(p,uX_scalar)),panel_ids)
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())
writevtk(Ω,dir*"/test_u",cellfields=["u"=>ucf],append=false)


# _u = map(p->GenericField(u_scalar(p)),panel_ids)
# ucf = CellField(u_scalar,Ω)

_
function met(p)

  function m(x)
    function y(x)
      if p == 1
        return x
      elseif p == 2
        return Ps[2] ⋅ (x- VectorValue(0,π/2))
      end
    end

    metric1((x))

  end
end

function invmet(p)
  function m(x)
    function y(x)
      if p == 1
        return x
      elseif p == 2
        return Ps[2] ⋅ (x- VectorValue(0,π/2))
      end
    end

    invmetric1((x))

  end
end

function sqrtmet(p)


  function m(x)
    function y(x)
      if p == 1
        return x
      elseif p == 2
        return Ps[2] ⋅ (x- VectorValue(0,π/2))
      end
    end

    sqrtmet1((x))

  end
end

pt = Point(π/4,π/4)



pt2 = Point(π/8,3π/8)
# pt2 = Point(π/4,3π/4)
pt1 = Ps[2] ⋅ (pt2 - VectorValue(0,π/2))

metric1(pt1)
metric2(pt2)
metric1(pt2)

invmetric1(pt1)
invmetric2(pt2)

sqrtmet1(pt1)
sqrtmet2(pt2)

_m = map(p->GenericField(met(p)),panel_ids)
M = CellData.GenericCellField(_m,Ω,PhysicalDomain())

# sqrtM(get_cell_points(Ω))

_minv = map(p->GenericField(invmet(p)),panel_ids)
invM = CellData.GenericCellField(_minv,Ω,PhysicalDomain())

_msqrt = map(p->GenericField(sqrtmet(p)),panel_ids)
sqrtM = CellData.GenericCellField(_msqrt,Ω,PhysicalDomain())

mNew =  Metric(M,sqrtM,invM,met,sqrtmet,invmet)

lapucf = 1/sqrtM * divergence( sqrtM*( invM ⋅ gradient(ucf) ) )
# (surface_laplacian(ucf,mNew)-lapucf)(pts)
rhs = ucf + lapucf
sum(∫(rhs)dΩ)
# mass(u,v) = ∫( (u*v)*sqrtmeasG )dΩ
# stiffnes(u,v) =∫( ( (Jcf ⋅(invG ⋅ ∇(v)))⋅ (Jcf ⋅invG ⋅ ∇(u) ) ) *(sqrtmeasG))dΩ
# bound(v) = ∫( (((invG ⋅∇(ucf))⋅nΓ*v)*sqrtmeasG ) )dΓ
# force(v) = ∫(  (rhs*v)*sqrtmeasG )dΩ

mass(u,v) = ∫( (u*v)*sqrtM )dΩ
stiffnes(u,v) =∫( (∇(v)⋅ (invM⋅ ∇(u) )) *sqrtM)dΩ
bound(v) = ∫( (((invM ⋅∇(ucf))⋅(nΓ)*v)*sqrtM ) )dΓ
force(v) = ∫(  (rhs*v)*sqrtM )dΩ

helmholtz_biform(u,v) = mass(u,v) - stiffnes(u,v)
helmholtz_liform(v) = force(v) - bound(v)


op = AffineFEOperator(helmholtz_biform,helmholtz_liform,U,V)

uh = solve(LUSolver(),op)

function _u_scalar_ambient2parametric2(p::Int,uX::Function)
  function _u(αβ)
    αβ1 =  ShiftField(p)(αβ)
    # _αβ1 = InversionField(Ps[p])(αβ1)
    θϕ = GnomonicField()(αβ1)
    _XYZ = SigmaField(RADIUS)(θϕ)
    XYZ = PanelRotationField(r1p_3D[p])(_XYZ)

    uX(XYZ)
  end
end

_cell_field = map(p->GenericField(_u_scalar_ambient2parametric2(p,uX_scalar)),panel_ids)
_ucf = CellData.GenericCellField(_cell_field,Ω,PhysicalDomain())

e = l2(uh-ucf,dΩ)

writevtk(Ω,dir*"/test_u",cellfields=["u"=>ucf,"uh"=>uh,"e"=>ucf-uh,"_u"=>_ucf,"_e"=>_ucf-uh,"eu"=>_ucf-ucf],append=false)

# writevtk(Γ,dir*"/flux12",cellfields=["n"=>nΓ,"nrot"=>Pcf⋅nΓ],append=false)
