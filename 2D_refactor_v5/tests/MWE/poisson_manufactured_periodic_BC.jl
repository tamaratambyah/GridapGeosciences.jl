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

using Gridap
# using DrWatson
# dir = datadir("2D_CubedSphereRefactor")
# !isdir(dir) && mkdir(dir)


#### Analytic solution
p = 2
u(x) = x[1]*(1-x[1])

# p = 4
# u(x) = x[1]*(1-x[1])*x[2]*(1-x[2])


#### Domain and quadrature
degree = 2*(p+1)

model = CartesianDiscreteModel((0,1, 0,1), (10,10), isperiodic=(true,true))
# writevtk(model,dir*"/model",append=false)

Ω = Triangulation(model)
dΩ = Measure(Ω,degree)


#### FE Problem
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

ucf = CellField(u,Ω)
u_analytic = interpolate(u,U)
# writevtk(Ω,dir*"/poisson_manufactured_periodic_BC",
#         cellfields=["u"=>u,"ucf"=>ucf,"uh"=>u_analytic ],append=false)

rhs = -1.0*laplacian(u_analytic)
poisson_biform(u,v) = ∫( gradient(u)⋅gradient(v)  )dΩ
poisson_liform(v) = ∫(  rhs*v )dΩ

op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

#### Compute errors
e = uh-u_analytic
sum(∫(e⊙e)dΩ)


# writevtk(Ω,dir*"/poisson_manufactured_periodic_BC",
#         cellfields=["u"=>u,"ucf"=>ucf,"uh"=>uh],append=false)
