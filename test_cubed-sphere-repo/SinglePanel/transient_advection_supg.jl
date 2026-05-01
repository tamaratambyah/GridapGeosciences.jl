using DrWatson
using Gridap
using GridapSolvers
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using GridapGeosciences
using Test
using DataFrames



include("../convergence_tools.jl")
include("../output_tools.jl")


u0(x) = sin(4*x[1])
velocity = VectorValue(-1,0)

model0 = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(4,4),isperiodic=(true,true))
model1 = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(8,8),isperiodic=(true,true))
model2 = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(16,16),isperiodic=(true,true))
model3 = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(32,32),isperiodic=(true,true))
models = [model3,model2,model1,model0]

_nc(model) = num_cells(model)


################################################################################
#### Transient
################################################################################
function transient_advection_supg_solver(
  panel_model,
  p_fe::Int,_dir::String,u::Function,v::Function,CFL=0.1,ls=LUSolver(),tF=2*π,return_vtk=false)


  lvl = nref(_nc(panel_model))
  println("nlevl = $lvl")

  dir = _dir*"/sol_p$(p_fe)_nref$lvl"
  ( !isdir(dir) ) && mkdir(dir)

  ## now enter the solver
  panel_ids = [1 for i in 1:num_cells(panel_model)]
  degree = 2*(p_fe + 1)

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)


  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  P = TransientTrialFESpace(Q)


  # supg stabilisation parameter
  _dx = (π/2)  / sqrt(num_cells(panel_model) )
  _dt = _dx*CFL#/p_fe^2
  # dt = floor(_dt,sigdigits=1)
  dt = _dt


  meas_cf = ParametricCellField(sqrtg,Ω_panel,panel_ids)


  β = CellField(velocity,Ω_panel)# VectorValue(-1,0 )
  τ = 0.5*dt

  a_Ω(t, u, v) = ∫( (v * β ⋅ ∇(u))*meas_cf )dΩ
  a_s(t, u, v) = ∫(( (β ⋅ ∇(u)) * (β ⋅ ∇(v)))*meas_cf)dΩ
  m_Ω(t, u, v) = ∫( (∂t(u) * v)*meas_cf)dΩ
  m_s(t, u, v) = ∫( (∂t(u) *(β ⋅ ∇(v)))*meas_cf)dΩ

  m(t, u, v) = m_Ω(t, u, v) + τ * m_s(t, u, v)
  a(t, u, v) = a_Ω(t, u, v) + τ * a_s(t, u, v)

  res(t,u,v) = m(t,u,v) + a(t,u,v)
  jac(t,u,du,v) = a(t,du,v)
  jac_t(t,u,dtu,v) = ∫( v*dtu )dΩ + τ * ∫( dtu *(β ⋅ ∇(v)))dΩ
  opT = TransientFEOperator(res, (jac, jac_t), P, Q)

  # solve with SSP RK 3
  uh0 = interpolate_everywhere(u0, P(0.0))
  t0 = 0.0

  nls = NLSolver(ls, show_trace=true, method=:newton, iterations=20)
  # solver = RungeKutta(nls, ls, dt, :EXRK_SSP_3_3)
  solver = RungeKutta(nls, ls, dt, :SDIRK_3_2)
  # solver = ThetaMethod(nls,dt,0.5)
  solT = solve(solver, opT, t0, tF, uh0)

  labels = ["uh", "eu"]

  Ω_error = Triangulation(panel_model)
  dΩ_error = Measure(Ω_error,6*p_fe+1)
  if return_vtk
    panel_cfs = [uh0, uh0-uh0]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=cellfields,append=false)
  end

  ## store errors
  ts = Float64[]
  Es = Float64[]
  Ms = Float64[]
  eu = 0.0
  t = 0.0

  push!(ts,0.0)
  push!(Es,0.0)

  counter = 1

  _uh = uh0
  for (t,uh) in solT

    println("t = ", t)

    eu = l2((uh-uh0),dΩ_error)

    push!(ts,t)
    push!(Es,eu)

    _uh = uh
    if return_vtk && (mod(counter,10) == 0)
      panel_cfs = [uh, uh-uh0]
      cellfields = map((x,y) -> x=>y, labels,panel_cfs)

      writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=cellfields,append=false)
    end
    counter = counter + 1
  end

  push!(ts,t)
  push!(Es,eu)

  panel_cfs = [_uh, _uh-uh0]
  cellfields = map((x,y) -> x=>y, labels,panel_cfs)
  writevtk(Ω_panel,_dir*"/final_sol_p$(p_fe)_nref$lvl.vtu", cellfields=cellfields,append=false)

  make_pvd(dir,"solT",1)

  ### convergence output for DrWatson
  dir_convergence = _dir*"/convergence"
  (!isdir(dir_convergence)) && mkdir(dir_convergence)

  n = nc(panel_model)
  dxx = dx(panel_model)
  output = @strdict n dxx p_fe lvl ts Es
  safesave(datadir(dir_convergence, ("transient_advection_dg_nref$(lvl)_p$p_fe.jld2")), output)


  return ts, Es
end

## helper function to return the error for transient solution
function transient_advection_supg_errors(panel_model,args...)
  ts, Es  = transient_advection_supg_solver(panel_model,args...)
  # return minimum(Es[end-10:end]),false,false
  return minimum(Es[end]),false,false
end



################################################################################
#### Convergence test with plots
################################################################################
ps = [1]
ls = LUSolver()
CFL = 0.5

v = u0
u = u0
tF = π/2

# u = panel_to_cartesian(u0)
# tF = 2*π

return_vtk = true

simName = "transient_advection_supg_convergence"

dir = datadir("TransientAdvectionSUPG_SinglePanel_sdirk")
(!isdir(dir)) && mkdir(dir)

errors = Vector{Vector{Float64}}(undef,length(ps))
ns = Vector{Vector{Float64}}(undef,length(ps))
dxs = Vector{Vector{Float64}}(undef,length(ps))
slopes = Vector{Float64}(undef,length(ps))

i = 1
p_fe = ps[i]
errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,transient_advection_supg_errors,p_fe,dir,u,v,CFL,ls,tF,return_vtk)


print_convergence_results(errors,ns,dxs,slopes,ps)

output = @strdict errors ns dxs slopes ps
safesave(datadir(dir, ("$simName.jld2")), output)

plot_convergence_from_saved(dir,simName)


##########
dir = datadir("TransientAdvectionSUPG_SinglePanel_sdirk/convergence")
df = collect_results(dir)
varNames = ["u"]
ps = unique(df[!,:p_fe])

plot()
for (i,p_fe) in enumerate(ps)
  E = df[(df.p_fe .== p_fe ),:Es]
  errors = sqrt.(map(x->minimum(x[end]), E))
  ns = df[(df.p_fe .== p_fe ),:n]
  dxs = (π/2)./ns
  slope = convergence_rate(dxs,errors)
  leginf = map(x->"$x: p=$p_fe", varNames)
  cols = map(x->palette(:tab10)[p_fe], varNames )
  plot_convergence(errors,ns,dxs,slope;leginf=leginf,colors=cols)
end
plot!(show=true)
savefig(dir*"/convergence_explicit")




E = df[(df.p_fe .== 3 ),:Es]
e3 = sqrt.(map(x->minimum(x[end-150:end]), E))


_linestyle = [:solid, :dash, :dot, :dashdotdot, :dashdot, :solid, :dash,]
_markers= [:circle, :rect, :star5, :diamond,  :cross, :hexagon]
_colors = palette(:tab10)

markers= [:circle :rect  :diamond ]
markersize = [6 7 6]
linestyle = [:solid :dash]
colors = [_colors[1] _colors[2] _colors[3]]

default(; fontfamily="Computer Modern");
plot(ns,[e2,e3],
    lw=2,
    marker=markers[1],
    ms = markersize[1],
    color=colors,
    label=[latexstring("\$ p = 2 \$") latexstring("\$ p = 3 \$")])
plot!(shape=:auto,
    xaxis=:log2,yaxis=:log10,
    xlabel=latexstring("\$ {n} \$"),
    ylabel=L"$L^2$ error",
    xtickfontsize=11,ytickfontsize=11,
    xguidefontsize=12,yguidefontsize=12,
    legendfontsize=10,
    legend=:bottomleft,
    #legend_columns=2,
    framestyle = :box,
    # guidefontfamily=font(20,"Times Roman")
    ylimits=(1e-6,0)
    )
