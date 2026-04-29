using Gridap
using GridapGeosciences
using DrWatson

include("Laplace/analytic_funcs.jl")
include("convergence_tools.jl")
include("output_tools.jl")

s_models = get_refined_models(3)

model = s_models[2]
Ω_panel = Triangulation(model)
panel_ids = get_panel_ids(model)
u_cf = ParametricCellField(f_sin,Ω_panel,panel_ids)


Q = TestFESpace(model, ReferenceFE(lagrangian,Float64,1); conformity=:L2)
P = TrialFESpace(Q)

uh = interpolate(u_cf,P)


save_mesh(plotsdir(),model)

t = 0.0
cellfields=[uh,uh]
labels = ["uh1", "uh2"]

save_cellfields(plotsdir(),Triangulation(model),t,cellfields,labels)

get_cell_dof_values(u_cf)



using ScatteredInterpolation
using GLMakie

include("plot_tools.jl")

dir = datadir("TransientAdvectionSUPG/sol_p1_nref4")
n = 100
plotName = "Advection"
cf = "uh"
plot_latlon(dir,n,plotName,cf)


######################### Distributed Model

using MPI
using PartitionedArrays
using GridapDistributed
using MPIPreferences
MPIPreferences.use_jll_binary()

nprocs = 6

ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end

include("convergence_tools.jl")
include("output_tools.jl")
models = get_distributed_refined_models(ranks,nprocs,3)

dir = plotsdir()
model = models[3]
save_mesh(dir,model)


Ω_panel = Triangulation(model)
panel_ids = get_panel_ids(model)
u_cf = ParametricCellField(f_sin,Ω_panel,panel_ids)

Q = TestFESpace(model, ReferenceFE(lagrangian,Float64,1); conformity=:L2)
P = TrialFESpace(Q)

uh = interpolate(u_cf,P)

cellfields = [uh, uh]
labels = ["uh", "ph"]
t = 100.0
save_cellfields(dir,Ω_panel,t,cellfields,labels)
