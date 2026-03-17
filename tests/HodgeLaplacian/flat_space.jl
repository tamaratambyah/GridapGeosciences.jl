using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using Gridap.ReferenceFEs, Gridap.Polynomials, Gridap.CellData

using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Test
using GridapPETSc

include("../convergence_tools.jl")

dir = datadir("HodgeLaplacian_flat")
!isdir(dir) && mkdir(dir)

model0 = CartesianDiscreteModel((0, 1, 0, 1, 0, 1), (2, 2, 2))
model1 = CartesianDiscreteModel((0, 1, 0, 1, 0, 1), (4, 4, 4))
model2 = CartesianDiscreteModel((0, 1, 0, 1, 0, 1), (16, 16, 16))
model3 = CartesianDiscreteModel((0, 1, 0, 1, 0, 1), (32, 32, 32))

models = [model3, model2, model1, model0]


function uX(x)
  VectorValue(x[1]^2,x[2]^2,x[3]^2)
end

divuX(x) = (∇⋅uX)(x)
ccurl(x) = (∇×(∇×uX))(x)
rhs(x) = ccurl(x) - gradient(divuX)(x)

u(x) = x[1]^2
f(x) = - Δ(u)(x)

panel_model = models[2]
writevtk(Triangulation(panel_model),dir*"/rhs",
           cellfields=["u"=>uX, "div"=>divuX, "rhs"=>rhs, "r1"=>f],append=false)

ls = LUSolver()
p_fe = 1


degree = 4*(p_fe + 1)
if p_fe == 0
  degree = 8
end

## finite element solver
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,2*degree)
Ω_error = Triangulation(panel_model)
dΩ_error = Measure(Ω_error,4*degree)

## cellfields
u_cf = CellField(uX,Ω_panel)
div_cf = CellField(divuX,Ω_panel)
rhs_cf = CellField(rhs,Ω_panel)
sigma_cf = -div_cf


tags = ["boundary"]

## FE spaces
T = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1,dirichlet_tags=tags)
S = TrialFESpace(T,sigma_cf)

R = TestFESpace(Ω_panel, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
H = TrialFESpace(R,u_cf)

X = MultiFieldFESpace([S,H])
Y = MultiFieldFESpace([T,R])

biform_s((s,u),(t,v)) = ∫( s*t )dΩ         - ∫( u⋅∇(t) )dΩ
biform_u((s,u),(t,v)) = ∫( (∇×u)⋅(∇×v) )dΩ + ∫( v⋅∇(s) )dΩ
biformX((s,u),(t,v)) = biform_s((s,u),(t,v)) + biform_u((s,u),(t,v))
liformX((t,v)) = ∫( rhs⋅v )dΩ
op = AffineFEOperator(biformX,liformX,X,Y)
sh, uh = solve(ls,op)

u_int = interpolate(u_cf,H)

_e =  uh - u_int
el2_u =  sqrt(sum(∫( _e⋅_e)dΩ_error))

_e = sh - sigma_cf
el2_s =  sqrt(sum(∫(  _e*_e )dΩ_error))


cellfields = ["u"=>u_int,
              "uh"=>uh,
              "eu"=>u_int-uh,
              "s"=>sigma_cf,
              "sh"=>sh,
              "es"=>sigma_cf-sh  ]

### plot in 3D
writevtk(Ω_panel,dir*"/flat_sol",
        cellfields=cellfields,append=false)
