using Gridap
using DrWatson
include("../convergence_tools.jl")

function H1projection(model,
  p_fe::Int,dir::String,f::Function,return_vtk)
  # 1. Define the geometry (a square domain)


  # 2. Define the Finite Element Space (H1-conforming, Lagrangian)
  reffe = ReferenceFE(lagrangian, Float64, p_fe)
  Vh = TestFESpace(model, reffe, conformity=:H1)
  Uh = TrialFESpace(Vh)

  # 3. Define the function we want to project
  # f(x) = sin(2*π*x[1]) * cos(2*π*x[2])
  # f(x) = sin(π*x[1]) * sin(π*x[2])
  # f(x) = x[2]*(1-x[2])
  # writevtk(Ω, "h1_projection", cellfields=[ "f"=>f],append=false)


  # 4. Define the Integration (Lebesgue measure)
  order = 1
  degree = 2*order+1
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)

  # 5. Define the H1 Inner Product (The Projection "Problem")
  # Bilinear form: <u, v>_H1
  a(u, v) = ∫( u*v + ∇(u)⋅∇(v) )*dΩ

  # Linear form: <f, v>_H1
  # l(v) = ∫( f*v + ∇(f)⋅∇(v) )*Measure(Ω, 15)
  l(v) = ∫( f*v -Δ(f)*v )*Measure(Ω, 15)
  sum(∫(f)dΩ)


  # 6. Solve for the projected function u_h
  op = AffineFEOperator(a, l, Uh, Vh)
  uh = solve(op)

  e=f-uh
  el2 = sum(∫( e*e )*dΩ)
  eh1 = sum(∫( e*e + ∇(e)⋅∇(e) )*dΩ)


  # 7. Visualize (Outputs a .vtu file for Paraview)
  writevtk(Ω, dir*"/h1_projection", cellfields=["uh"=>uh, "error"=>f-uh, "f"=>f],append=false)

  return eh1, false, false
  # return el2, false, false
end


model1 = CartesianDiscreteModel((0, 1, 0, 1), (4, 4),isperiodic=(true,true))
model2 = CartesianDiscreteModel((0, 1, 0, 1), (16, 16),isperiodic=(true,true))
model3 = CartesianDiscreteModel((0, 1, 0, 1), (32, 32),isperiodic=(true,true))
model4 = CartesianDiscreteModel((0, 1, 0, 1), (64, 64),isperiodic=(true,true))

models = [model4, model3, model2, model1]

# f(x) = sin(2*π*x[1]) * cos(2*π*x[2])
# f(x) = sin(π*x[1]) * sin(π*x[2])
f(x) = x[2]*(1-x[2])


ps = [1]

dir = datadir("InterpolationConvergence")
!isdir(dir) && mkdir(dir)

_dir = dir*"/scalar_func_2D_H1proj"
!isdir(_dir) && mkdir(_dir)
p_convergence_test([true],ps,models,H1projection,_dir,f,true)
plot_convergence_from_saved(_dir,"convergence",["H1 norm","L2 norm", ])
