"""
Manufacture solutions for Possion problem on a flat domain with periodic BC.
This test is a precursor to manufacturing solutions on the sphere.

Solve Δu = -f on Ω = [0,1]^2 doubly periodic.
The weak form is:
  ∫ (∇u)⋅(∇v) dΩ = ∫ fv dΩ
where f = -Δ(u_ex) is the rhs data.

Since the domain is doubly periodic, there are no boundary integrals in the weak
form. In addition, consider u_ex as periodic functions on [0,1]^2.
Examples are:
  u_ex = x(1-x) for p = 2 polynomial
      * this function is zero on x=0, x=1, and non-zero (periodic) on y=0, y=1

  u_ex = x(1-x)*y(1-y) for p = 4 polynomial
      * this function is zero on x=0, x=1, y=0, y=1

Find that the error blows up to 1e28
"""

using Gridap, Gridap.Geometry
using LinearAlgebra
using DrWatson
dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


#### Analytic solution2
p = 2
u(x) = x[1]*(1-x[1])

# p = 4
# u(x) = x[1]*(1-x[1])*x[2]*(1-x[2])

_u(x) = u(x) - sum(∫( u )dΩ)


#### Domain and quadrature
degree = 2*(p+1)

model = CartesianDiscreteModel((0,1, 0,1), (3,3), isperiodic=(true,true))


# _model = CartesianDiscreteModel((0,1, 0,1), (10,10))

writevtk(model,dir*"/model",append=false)

Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

## assess zero mean of analytic function
sum(∫(u)dΩ)
sum(∫(_u)dΩ)


#### FE Problem -- no lagrange multiplers
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,constraint=:zeromean)
U = TrialFESpace(V)

ucf = CellField(_u,Ω)
u_analytic = interpolate(_u,U)
sum(∫(u_analytic)dΩ)

# writevtk(Triangulation(_model),dir*"/poisson_manufactured_periodic_BC",
#         cellfields=["u"=>u ],append=false)


_rhs(x) = -1.0*(laplacian(_u)(x))
rhs(x) = -1.0*(laplacian(u)(x))

# poisson_biform(u,v) = ∫( gradient(u)⋅gradient(v)  )dΩ
# poisson_liform(v) = ∫(  _rhs*v )dΩ
# op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
# uh = solve(LUSolver(),op)

#### FE Problem -- with lagrange multiplers to impose zero mean
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)
Λ = ConstantFESpace(model)
M = TrialFESpace(Λ)
X = MultiFieldFESpace([U,M])
Y = MultiFieldFESpace([V,Λ])

poisson_biform((u,μ),(v,λ)) = ∫( gradient(u)⋅gradient(v)  )dΩ  + ∫(v*μ)dΩ + ∫(λ*u)dΩ
poisson_liform((v,λ)) = ∫(  _rhs*v )dΩ + ∫(λ*_u)dΩ

op = AffineFEOperator(poisson_biform,poisson_liform,X,Y)
uh,μh = solve(LUSolver(),op)

sum( ∫(uh )dΩ  )
sum( ∫(_u )dΩ  )

# A = get_matrix(op)
# println(eigvals(Array(A)))
# cond(Matrix(A))

#### Compute errors
e = uh-_u
sum(∫(e⊙e)dΩ)

# res((v,λ)) = poisson_biform((u,2/3),(v,λ))  - poisson_liform((v,λ))
# _r = assemble_vector(res,Y)
# norm(_r)



res2(v) = ∫( gradient(u)⋅gradient(v)  )dΩ  -  ∫(  _rhs*v )dΩ
_r = assemble_vector(res2,V)
norm(_r)



sum( ∫(_u )dΩ )

writevtk(Ω,dir*"/poisson_manufactured_periodic_BC",
        cellfields=["u"=>u,"uu"=>_u, "uh"=>uh, "e"=>uh-_u],append=false)


findall(x -> abs(x) > 1e-8,_r)
get_cell_dof_ids(V)
