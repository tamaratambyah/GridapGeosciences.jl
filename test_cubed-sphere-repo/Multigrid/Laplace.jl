using DrWatson
using PartitionedArrays
using Gridap, Gridap.Algebra
using Gridap.Helpers, Gridap.Adaptivity
using Test

using GridapSolvers
using GridapSolvers.MultilevelTools
using GridapSolvers.LinearSolvers
using GridapGeosciences

# using MPI, PartitionedArrays, GridapP4est, GridapDistributed
using Plots

include("../Laplace/analytic_funcs.jl")
include("../convergence_tools.jl")

models = get_refined_models(3,true)
mh = ModelHierarchy(models)

# FE Spaces
p_fe  = 1
qdegree = 2*p_fe + 1
reffe  = ReferenceFE(lagrangian,Float64,p_fe)
tests  = TestFESpace(mh,reffe;conformity=:H1)
trials = TrialFESpace(tests)


setup_prolongation_operators(trials,qdegree;)




f = panel_to_cartesian(fθϕ)
