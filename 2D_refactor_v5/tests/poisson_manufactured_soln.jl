using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots, LaTeXStrings
include("../src/initialise.jl")

_colors = palette(:tab10)
markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross]
ls = [:solid, :dash, :dot, :dashdot, :dashdotdot]

function solve_poisson(domain,n,p,degree,u,metric_func;isperiodic=fill(false,length(n)))

  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,n,isperiodic=isperiodic))

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

function plot_error(ns,errs,leginf,colors,ls,markers)
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

function plot_convergence(ns,errs,errs_g,leginf)
  plot()
  plot_error(ns,errs,leginf,_colors,
    fill(:solid,length(leginf)),fill(:circle,length(leginf)))

  plot_error(ns,errs_g,fill(false,length(leginf)),_colors,
  fill(:dash,length(leginf)),fill(:xcross,length(leginf)))

  plot!(yscale=:log10,framestyle=:box,
  xscale=:log10,
  xlabel="n cells",
  ylabel=L"L2(u_{ex} - u_h)"
  )
  plot!(show=true)
  plot!(xtickfontsize=10,ytickfontsize=10,
  legendfontsize=10,guidefontsize=10)
end


p = 2
degree = 2*(p+1)
ns = [2^i for i = 1:6]
dx = 1 ./ ns

################################################################################
#### 1D tests
################################################################################
domain = (0,1)
# u(x) = x[1]*(1-x[1])
u(x) = cos(x[1])
uambient(x) = cos(x[1])*sin(x[2]) + 2*x[1]^2

# dictionary of metric, uambient(geomap)
dd = Dict("linear" => ( l(r) = x -> if (r) TensorValue{1}(4)
                                    else uambient(VectorValue(x[1],2*x[1])) end ),
          "quad" => ( q(r) = x -> if (r) TensorValue{1}(4*x[1]^2 + 1)
                                  else uambient(VectorValue(x[1],x[1]^2)) end ),
          "cubic" => ( c(r) = x -> if (r) TensorValue{1}(9*x[1]^4 + 12*x[1]^3 + 10*x[1]^2 + 4*x[1] + 2)
                                  else uambient( VectorValue(x[1],x[1]^3+x[1]^2+x[1]+1)) end )
 )

# convergence plot - exact solution in parametric space
errs = []
errs_g = []
for (key, val) in dd
  for n in collect(ns)
    e, eg = solve_poisson(domain,(n),p,degree,u,val(true))
    push!(errs,e)
    push!(errs_g,eg)
  end
end
plot_convergence(ns,errs,errs_g,collect(keys(dd)))
# savefig(plotsdir()*"/poisson_parametric_manufactured_1D")
savefig(plotsdir()*"/poisson_parametric_convergence_1D")

# convergence plot - exact solution in ambient space
errs = []
errs_g = []
for (key, val) in dd
  for n in collect(ns)
    e, eg = solve_poisson(domain,(n),p,degree,(val(false)),val(true))
    push!(errs,e)
    push!(errs_g,eg)
  end
end
plot_convergence(ns,errs,errs_g,collect(keys(dd)))
savefig(plotsdir()*"/poisson_ambient_convergence_1D")


################################################################################
#### 2D tests
################################################################################
domain = (0,1,0,1)
u(x) = x[1]*(1-x[1]) + x[2]*(1-x[2])
# u(x) = cos(x[1])*sin(x[2])
uambient(x) = cos(x[1])*sin(x[2])*x[3]


dd = Dict(
# "const" => ( cc(r) = x -> if (r) TensorValue{2,2}(1,0,0,1 )
                                    # else uambient(VectorValue(x[1],x[2],0.0)) end ),
          "linear" => ( ll(r) = x -> if (r) TensorValue{2,2}(2,1,1,2 )
                                    else uambient(VectorValue(x[1],x[2],x[1]+x[2])) end ),
          "quad" => ( qq(r) = x -> if (r) TensorValue{2,2}(1+4*x[1]^2,4*x[1]*x[2],4*x[1]*x[2],1+4*x[2]^2 )
                                  else uambient(VectorValue(x[1],x[2],x[1]^2+x[2]^2)) end )
 )


# convergence plot - exact solution in parametric space
errs = []
errs_g = []
for (key, val) in dd
  for n in collect(ns)
    e, eg = solve_poisson(domain,(n,n),p,degree,u,val(true))
    push!(errs,e)
    push!(errs_g,eg)
  end
end
plot_convergence(ns,errs,errs_g,collect(keys(dd)))
savefig(plotsdir()*"/poisson_parametric_manufactured_2D")
# savefig(plotsdir()*"/poisson_parametric_convergence_2D")



# convergence plot - exact solution in ambient space
errs = []
errs_g = []
for (key, val) in dd
  for n in collect(ns)
    e, eg = solve_poisson(domain,(n,n),p,degree,(val(false)),val(true))
    push!(errs,e)
    push!(errs_g,eg)
  end
end
plot_convergence(ns,errs,errs_g,collect(keys(dd)))
savefig(plotsdir()*"/poisson_ambient_convergence_2D")

################################################################################
#### Cylinder tests
################################################################################
ns = [2^i for i = 2:6]

u(x) = x[2]*(1-x[2])
# u(x) = cos(x[1])
domain = (0,2*π,0,1.0)
metric_func(x) = TensorValue{2,2}(1,0.0,0.0,1.0)

errs = []
errs_g = []
for n in collect(ns)
  e, eg = solve_poisson(domain,(n,n),p,degree,u,metric_func;isperiodic=(true,false))
  push!(errs,e)
  push!(errs_g,eg)
end

plot_convergence(ns,errs,errs_g,["cylinder"])

dx = (2π ./ ns ) .* (1 ./ ns)
gg = [5e-3dx.^3]
plot!(ns,gg[1],lw=2,c=:black,ls=:dash,label=latexstring("\$ (\\Delta x)^3 \$"))

savefig(plotsdir()*"/poisson_parametric_manufactured_2D_cylinder")
# savefig(plotsdir()*"/poisson_parametric_convergence_2D_cylinder")





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
### Sphere -- not working

p = 4
degree = 2*(p+1)
u(x) = x[1]*(1-x[1])*x[2]*(1-x[2]) #(x[1]-π/2)*(x[1]+π/2) + (x[2]-2*π)*x[2]
# u(x) = cos(x[1])

# model = CartesianDiscreteModel((-π/2,π/2, 0,2*π ),(10,10),isperiodic=(true,true))
model = CartesianDiscreteModel((0,1,0,1 ),(10,10),isperiodic=(true,true))

writevtk(model,dir*"/sphere_model",append=false)

Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

metric_func(x) = TensorValue{2,2}(1,0.0,0.0,(cos(x[1]))^2 )

m = Metric(metric_func,Ω)
dΩg =  Measure(m,Ω,degree)

#### FE Problem

V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)


ucf = CellField(u,Ω)
u_analytic = interpolate(u,U)
writevtk(Ω,dir*"/sphere",
        cellfields=["u"=>u,"ucf"=>ucf,"uh"=>u_analytic ],append=false)

# rhs = -1.0*surface_laplacian(ucf,m)
# poisson_biform(u,v) = ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
# poisson_liform(v) = ∫(  rhs*v )dΩg

rhs = -1.0*laplacian(ucf)
poisson_biform(u,v) = ∫( gradient(u)⋅gradient(v)  )dΩ
poisson_liform(v) = ∫(  rhs*v )dΩ

op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(uh-u_analytic,dΩ)
l2(uh-u_analytic ,dΩg)
writevtk(Ω,dir*"/sphere",
        cellfields=["u"=>u,"ucf"=>ucf,"uh"=>uh],append=false)
