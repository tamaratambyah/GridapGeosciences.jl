using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
include("../src/initialise.jl")


function solve_poisson(p,degree,u,model,metric_func)

  Ω = Triangulation(model)
  m = Metric(metric_func,Ω)

  dΩ = Measure(Ω,degree)
  dΩg =  Measure(m,Ω,degree)

  #### FE Problem

  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,
                        dirichlet_tags="boundary")
  U = TrialFESpace(V,u)

  ucf = CellField(u,Ω)

  rhs = -1.0*surface_laplacian(ucf,m)

  poisson_biform(u,v) = ∫( surface_gradient(u,m)⋅surface_gradient(v,m) )dΩg
  # poisson_liform(v) = ∫( (rhs*v) )dΩg
  poisson_liform(v) = ∫(  Operation(meas)(m.inv_metric)*rhs*v )dΩg

  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
  uh = solve(LUSolver(),op)

  println("Errors: ", l2(uh-ucf,dΩ), "; ", l2(uh-ucf,dΩg))
  return uh, ucf
end


p = 2
degree = 2*(p+1)

u(x) = x[1]*(1-x[1])
model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1),(1)))

metric_func(x) = TensorValue{1}(4)
solve_poisson(p,degree,u,model,metric_func)


########################### low level
### 1D tests

u(x) = x[1]*(1-x[1])

model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1),(1)))
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

metric_func(x) = TensorValue{1}(4)
# metric_func(x) = TensorValue{1}(4*x[1]^2 + 1)
# metric_func(x) = TensorValue{1}(9*x[1]^4 + 12*x[1]^3 + 10*x[1]^2 + 4*x[1] + 2)

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
poisson_biform(u,v) = ∫( surface_gradient(u,m)⋅surface_gradient(v,m) )dΩg
# poisson_liform(v) = ∫( Operation(meas)(m.inv_metric)*(rhs*v) )dΩg
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
# m = Metric(metric_func,Ω)

metric_func(x) = TensorValue{2,2}(2,1,1,2 )
m = Metric(metric_func,Ω)

# metric_func(x) = TensorValue{2,2}(1+4*x[1]^2,4*x[1]*x[2],4*x[1]*x[2],1+4*x[2]^2 )
# m = Metric(metric_func,Ω)

dΩg =  Measure(m,Ω,degree)

#### FE Problem

V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,
                      dirichlet_tags="boundary")
U = TrialFESpace(V,u)

ucf = CellField(u,Ω)

rhs = -1.0*surface_laplacian(ucf,m)

# poisson_biform(u,v) = ∫( ∇(u)⋅∇(v) )dΩ
# poisson_liform(v) = ∫( rhs*v )dΩ
poisson_biform(u,v) = ∫( surface_gradient(u,m)⋅surface_gradient(v,m) )dΩg
poisson_liform(v) = ∫(  Operation(meas)(m.inv_metric)*rhs*v )dΩg

op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(uh-ucf,dΩ)
l2(uh-ucf,dΩg)
