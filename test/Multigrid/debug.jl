using MPI, PartitionedArrays
using GridapGeosciences
using Gridap
using GridapP4est
using GridapDistributed


MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
ranks = distribute_with_mpi(LinearIndices((prod(nprocs),)))


o3model = Parametric3DOctreeDistributedDiscreteModel(ranks;
num_horizontal_uniform_refinements=0, num_vertical_uniform_refinements=0);

# level 1
o3model, = adapt_model(ranks,o3model)


mh = ModelHierarchy(o3model,2)
@test isa(mh,ModelHierarchy)

n_ref_lvls = 2
coarse_model = o3model
ranks = get_parts(coarse_model.parametric_dmodel)

models = Vector{GridapDistributed.DistributedDiscreteModel}(undef,n_ref_lvls+1)
# glues = Vector{MPIArray}(undef,n_ref_lvls)

models[n_ref_lvls+1] = coarse_model.parametric_dmodel

model = coarse_model
for n in n_ref_lvls:-1:1
  model, glue = adapt_model(ranks,model)
  models[n] = model.parametric_dmodel
  # glues[n] = glue
end

nlevs = length(models)
lev = 1
model1 = models[lev]
model2 = models[lev+1]
Gridap.Adaptivity.is_child(models[lev+1],models[lev])


 GridapSolvers.MultilevelTools.ModelHierarchy(models)
