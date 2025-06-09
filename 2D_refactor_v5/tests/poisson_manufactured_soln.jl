using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots, LaTeXStrings
include("../src/initialise.jl")

_colors = palette(:tab10)
markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross]
ls = [:solid, :dash, :dot, :dashdot, :dashdotdot]

function solve_poisson(domain,n,p,degree,u,metric_func)

  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,n))

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

  poisson_biform(u,v) = ∫( surface_gradient(u,m)⋅gradient(v) )dΩg
  poisson_liform(v) =  ∫( (rhs*v) )dΩg

  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
  uh = solve(LUSolver(),op)

  e = l2(uh-ucf,dΩ)
  eg = l2(uh-ucf,dΩg)

  println("Errors: ", e, "; ", eg)
  return e, eg
end

function plot_convergence(ns,errs,leginf,colors,ls,markers)
  nsims = Int(length(errs)/length(ns))
  for i in 1:nsims
    idx1 = 1 + (i-1)*length(ns)
    idx2 = (i)*length(ns)
    plot!(ns,
          errs[idx1:idx2],
          lw=3,
          c=colors[i],ls=ls[i], markershape=markers[i],
    label=leginf[i])
  end

end

p = 2
degree = 2*(p+1)


######## 1D tests
u(x) = x[1]*(1-x[1])
domain = (0,1)
ns = [2^i for i = 1:4]

metric_func1(x) = TensorValue{1}(4)
metric_func2(x) = TensorValue{1}(4*x[1]^2 + 1)
metric_func3(x) = TensorValue{1}(9*x[1]^4 + 12*x[1]^3 + 10*x[1]^2 + 4*x[1] + 2)

ms = [metric_func1, metric_func2, metric_func3]

errs = []
errs_g = []
for metric_func in ms
  for n in collect(ns)
    e, eg = solve_poisson(domain,(n),p,degree,u,metric_func)
    push!(errs,e)
    push!(errs_g,eg)
  end
end



leginf = ["Linear", "Quad", "Cubic"]

plot()
plot_convergence(ns,errs,leginf,_colors,
  fill(:solid,length(leginf)),fill(:circle,length(leginf)))
plot_convergence(ns,errs_g,
fill(false,length(leginf)),_colors,fill(:dash,length(leginf)),fill(:xcross,length(leginf)))
plot!(yscale=:log10,framestyle=:box,
# title = "surface area of cubed sphere",
xlabel="n cells",
ylabel=L"L2(u_{ex} - u_h)"
)
plot!(show=true)
plot!(xtickfontsize=10,ytickfontsize=10,
legendfontsize=10,guidefontsize=10)
savefig(plotsdir()*"/poisson_manufactured_1D")



######## 2D tests
u(x) = x[1]*(1-x[1]) + x[2]*(1-x[2])
domain = (0,1,0,1)
ns = [2^i for i = 1:4]

metric_func1(x) = TensorValue{2,2}(1,0,0,1 )
metric_func2(x) = TensorValue{2,2}(2,1,1,2 )
metric_func3(x) = TensorValue{2,2}(1+4*x[1]^2,4*x[1]*x[2],4*x[1]*x[2],1+4*x[2]^2 )


ms = [metric_func1, metric_func2, metric_func3]

errs = []
errs_g = []
for metric_func in ms
  for n in collect(ns)
    e, eg = solve_poisson(domain,(n,n),p,degree,u,metric_func)
    push!(errs,e)
    push!(errs_g,eg)
  end
end

leginf = ["Const", "Linear", "Quad"]
plot()
plot_convergence(ns,errs,leginf,_colors,
  fill(:solid,length(leginf)),fill(:circle,length(leginf)))
plot_convergence(ns,errs_g,
fill(false,length(leginf)),_colors,fill(:dash,length(leginf)),fill(:xcross,length(leginf)))
plot!(yscale=:log10,framestyle=:box,
# title = "surface area of cubed sphere",
xlabel="n cells",
ylabel=L"L2(u_{ex} - u_h)"
)
plot!(show=true)
plot!(xtickfontsize=10,ytickfontsize=10,
legendfontsize=10,guidefontsize=10)
savefig(plotsdir()*"/poisson_manufactured_2D")





########################### low level
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
metric_func(x) = TensorValue{2,2}(2,1,1,2 )
# metric_func(x) = TensorValue{2,2}(1+4*x[1]^2,4*x[1]*x[2],4*x[1]*x[2],1+4*x[2]^2 )


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
