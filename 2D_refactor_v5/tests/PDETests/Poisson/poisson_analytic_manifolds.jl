"""
Manufacture solutions for Possion problem on flat parametric spaces that is
associated to a manifold.

Solve Δu = -f on Ω, where
  - Δ is the Laplace-Beltrami operator,
  - Ω is the parametric space of manifold M
  - f = -surf_lap(u_ex) is a manufactured rhs

The weak form in terms of surface operators is:
  ∫ (surf_grad(u)) ⋅ (grad(v)) dΩg = ∫ v(-surf_lap(u_ex)) dΩg
where
  - ∫ is in the parametric space
  - surf_grad, surf_lap, dΩg account for the metric
  - grad is the standard flat gradient in the parametric space

Consider various mappings:
  1D interval -> polynomial: Ω = [0,1]; u_ex = x(1-x)
    * linear: φ(x) = (x,2x)
    * quad:   ϕ(x) = (x,x^2)
    * cubic:  ϕ(x) = (x,x^3 + x^2 + x + 1)

  2D square -> 3D plane: Ω = [0,1]^2; u_ex =  x(1-x) + y(1-y)
    * const:  φ(x,y) = (x,y,0)
    * linear: φ(x,y) = (x,y,x+y)
    * quad:   φ(x,y) = (x,y,x^2 + y^2)

  2D square -> cylinder: Ω = [0,2π] × [0,1]; u_ex = y(1-y)
    * 1 periodic boundary in parametric space

  2D square -> sphere: Ω = [-π/2, π/2 ] × [0,2π]
    * 2 periodic boundaries, requires zeromean constraint and function
"""



using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots, LaTeXStrings
include("../../../src/initialise.jl")
include("../pde_helpers.jl")
include("../analytic_metrics.jl")
include("PoissonSolvers.jl")

p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

################################################################################
#### 1D tests
################################################################################

# u(x) = x[1]*(1-x[1])
u(x) = cos(x[1])
uambient(x) = cos(x[1])*sin(x[2]) + 2*x[1]^2

# compute errors:
errs = []
errs_g = []

amb_errs = []
amb_errs_g = []

for name in manifolds_1D
  for n in collect(ns)

    # exact solution from parametric space
    e,eg = solve_poisson_manifold(name,n,p,degree,u)
    push!(errs,e)
    push!(errs_g,eg)

    # exact solution from ambient space
    amb_e, amb_eg = solve_poisson_manifold(name,n,p,degree,x->uambient(charts[name](x)))
    push!(amb_errs,amb_e)
    push!(amb_errs_g,amb_eg)
  end
end

leginf = map(x->string(x),collect(manifolds_1D))

# convergence plot - exact solution in parametric space
plot()
plot_error(ns,errs;leginf=leginf)
plot_error(ns,errs_g)
plot!(yscale=:log10,framestyle=:box,
  xscale=:log10,xlabel="n cells",ylabel=L"L2(u_{ex} - u_h)"
  )
plot!(ns,1e-4dx.^6,lw=2,c=:black,ls=:dash,label="dx^6")
savefig(plotsdir()*"/poisson_parametric_convergence_1D")

# convergence plot - exact solution in ambient space
plot()
plot_error(ns,amb_errs;leginf=leginf)
plot_error(ns,amb_errs_g)
plot!(yscale=:log10,framestyle=:box,
  xscale=:log10,xlabel="n cells",ylabel=L"L2(u_{ex} - u_h)"
  )
plot!(ns,1e-4dx.^6,lw=2,c=:black,ls=:dash,label="dx^6")
savefig(plotsdir()*"/poisson_ambient_convergence_1D")


################################################################################
#### 2D tests
################################################################################

# u(x) = x[1]*(1-x[1]) + x[2]*(1-x[2])
u(x) = cos(x[1])*sin(x[2])
uambient(x) = cos(x[1])*sin(x[2])*x[3]


# compute errors:
errs = []
errs_g = []
amb_errs = []
amb_errs_g = []
for name in [:linear2D,:quad2D]#manifolds_2D
  for n in collect(ns)

    # exact solution from parametric space
    e,eg = solve_poisson_manifold(name,n,p,degree,u)
    push!(errs,e)
    push!(errs_g,eg)

    # exact solution from ambient space
    amb_e, amb_eg = solve_poisson_manifold(name,n,p,degree,x->uambient(charts[name](x)))
    push!(amb_errs,amb_e)
    push!(amb_errs_g,amb_eg)
  end
end

leginf = map(x->string(x),collect([:linear2D,:quad2D]))

# convergence plot - exact solution in parametric space
plot()
plot_error(ns,errs;leginf=leginf)
plot_error(ns,errs_g)
plot!(yscale=:log10,framestyle=:box,
  xscale=:log10,xlabel="n cells",ylabel=L"L2(u_{ex} - u_h)"
)
# savefig(plotsdir()*"/poisson_parametric_manufactured_2D")
plot!(ns,1e-4dx.^6,lw=2,c=:black,ls=:dash,label="dx^6")
savefig(plotsdir()*"/poisson_parametric_convergence_2D")


# convergence plot - exact solution in ambient space
plot()
plot_error(ns,amb_errs;leginf=leginf)
plot_error(ns,amb_errs_g)
plot!(yscale=:log10,framestyle=:box,
  xscale=:log10,xlabel="n cells",ylabel=L"L2(u_{ex} - u_h)"
  )
plot!(ns,1e-4dx.^6,lw=2,c=:black,ls=:dash,label="dx^6")
savefig(plotsdir()*"/poisson_ambient_convergence_2D")

################################################################################
#### Cylinder tests
################################################################################
# u(x) = x[2]*(1-x[2])
u(x) = cos(x[1])

errs = []
errs_g = []
for n in collect(ns)
  e, eg =  solve_poisson_manifold(:cylinder,n,p,degree,u)
  push!(errs,e)
  push!(errs_g,eg)
end

plot()
plot_error(ns,errs)
plot_error(ns,errs_g)
plot!(yscale=:log10,framestyle=:box,
  xscale=:log10,xlabel="n cells",ylabel=L"L2(u_{ex} - u_h)"
)
# savefig(plotsdir()*"/poisson_parametric_manufactured_2D_cylinder")
dx = (2π ./ ns ) .* (1 ./ ns)
plot!(ns,5e-3dx.^3,lw=2,c=:black,ls=:dash,label=latexstring("\$ (\\Delta x)^3 \$"))
savefig(plotsdir()*"/poisson_parametric_convergence_2D_cylinder")



################################################################################
#### low level tests
################################################################################
###########################
### 1D tests

u(x) = x[1]*(1-x[1])

model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1),(10)))
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

# metric_func(x) = TensorValue{1}(4)
# metric_func(x) = TensorValue{1}(4*x[1]^2 + 1)
metric_func(x) = TensorValue{1}(9*x[1]^4 + 12*x[1]^3 + 10*x[1]^2 + 4*x[1] + 2)

m = Metric(metric_func,Ω)

dΩg =  Measure(m,Ω,degree)

#### FE Problem

V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,
                      dirichlet_tags="boundary")
U = TrialFESpace(V,u)

ucf = CellField(u,Ω)

rhs = -1.0*surface_laplacian(ucf,m)

# poisson_biform(u,v) = ∫( ∇(u)⋅∇(v) )dΩ
# poisson_liform(v) = ∫( rhs*v )dΩ
poisson_biform(u,v) = ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
poisson_liform(v) = ∫( (rhs*v) )dΩg

op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(uh-ucf,dΩ)
l2(uh-ucf,dΩg)



################################################################################
### 2D tests

u(x) = x[1]*(1-x[1]) + x[2]*(1-x[2])


model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(10,10)))
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

# metric_func(x) = TensorValue{2,2}(1,0,0,1 )
# metric_func(x) = TensorValue{2,2}(2,1,1,2 )
metric_func(x) = TensorValue{2,2}(1+4*x[1]^2,4*x[1]*x[2],4*x[1]*x[2],1+4*x[2]^2 )


m = Metric(metric_func,Ω)

dΩg =  Measure(m,Ω,degree)

#### FE Problem

V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,
                      dirichlet_tags="boundary")
U = TrialFESpace(V,u)


ucf = CellField(u,Ω)

rhs = -1.0*surface_laplacian(ucf,m)

# poisson_biform(u,v) = ∫( ∇(u)⋅∇(v) )dΩ
# poisson_liform(v) = ∫( rhs*v )dΩ
poisson_biform(u,v) = ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
poisson_liform(v) = ∫(  rhs*v )dΩg

op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(uh-ucf,dΩ)
l2(uh-ucf,dΩg)




################################################################################
### Cylinder

degree = 2
u(x) = x[2]*(1-x[2])
u(x) = cos(x[1])

model = CartesianDiscreteModel((0,2*π, 0,1),(10,10),isperiodic=(true,false))

writevtk(model,dir*"/cyclinder_model",append=false)

Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

metric_func(x) = TensorValue{2,2}(1,0.0,0.0,1.0 )

m = Metric(metric_func,Ω)
dΩg =  Measure(m,Ω,degree)

#### FE Problem

V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1, dirichlet_tags="boundary")
U = TrialFESpace(V,u)


ucf = CellField(u,Ω)

writevtk(Ω,dir*"/cyclinder",
        cellfields=["u"=>u,"ucf"=>ucf],append=false)

rhs = -1.0*surface_laplacian(ucf,m)
poisson_biform(u,v) = ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
poisson_liform(v) = ∫(  rhs*v )dΩg

# rhs = -1.0*laplacian(ucf)
# poisson_biform(u,v) = ∫( gradient(u)⋅gradient(v)  )dΩ
# poisson_liform(v) = ∫(  rhs*v )dΩ

op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(uh-ucf,dΩ)
l2(uh-ucf,dΩg)




################################################################################
### Sphere -- need to consider zero mean constraint in FE space and zero mean
### analytic solution

p = 2
degree = 2*(p+1)

function uex(x)
  if x[1] < 0
    return x[1]*(x[1]+π/2)
  else
    return  x[1]*(π/2-x[1])
  end
end

# uex(x) = cos(x[2] )
uex(x) = cos(2*x[1] )

model = CartesianDiscreteModel((-π/2,π/2, 0,2*π ),(10,10),isperiodic=(true,true))

# writevtk(model,dir*"/sphere_model",append=false)

Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
sum(∫( uex )dΩ ) # check zero mean

metric_func(x) = TensorValue{2,2}(1,0.0,0.0,(cos(x[1]))^2 )

m = Metric(metric_func,Ω)
dΩg =  Measure(m,Ω,degree)


#### FE Problem
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,constraint=:zeromean)
U = TrialFESpace(V)


ucf = CellField(uex,Ω)
sum( ∫( ucf )dΩ   )
writevtk(Ω,dir*"/sphere", cellfields=["u"=>uex,"ucf"=>ucf],append=false)

rhs = -1.0*surface_laplacian(ucf,m)
poisson_biform(u,v) = ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
poisson_liform(v) = ∫(  rhs*v )dΩg

op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(uh-ucf,dΩ)
l2(uh-ucf ,dΩg)
writevtk(Ω,dir*"/sphere",
        cellfields=["u"=>u,"ucf"=>ucf,"uh"=>uh],append=false)




### forcing zero mean
model = CartesianDiscreteModel((0,1,0,1 ),(4,4))
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)


uex(x) = x[1]*(1-x[1])
_uex(x) = uex(x) - sum(∫( uex )dΩ)

sum(∫( uex )dΩ ) # check zero mean
sum(∫( _uex )dΩ )

writevtk(Ω,dir*"/sphere",
        cellfields=["u"=>uex,"u0"=>_uex],append=false)
