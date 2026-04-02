using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using GridapGeosciences
using GridapP4est
using MPI
using PartitionedArrays
using Test

n_ref_lvls = 3
include("../convergence_tools.jl")

### serial test
models = get_refined_models(n_ref_lvls,true)

mh_serial = ModelHierarchy(models)
@test isa(mh_serial,ModelHierarchy)
reffe = ReferenceFE(lagrangian,Float64,1)
tests  = TestFESpace(mh_serial,reffe)


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

### distributed test
model0 = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=0)
mh = ModelHierarchy(model0,n_ref_lvls)
@test isa(mh,ModelHierarchy)
reffe = ReferenceFE(lagrangian,Float64,1)
tests  = TestFESpace(mh,reffe)
