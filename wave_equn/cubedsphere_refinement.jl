using GridapGeosciences
using GridapP4est
using PartitionedArrays
using MPI
using Gridap
using GridapDistributed
using DrWatson

function f(x)
    sin(2*pi*x[1])*cos(2*pi*x[2])*cos(2*pi*x[3])
end

function calculate_error_indicators(model,fh)
    Ω=Triangulation(model)
    dΩ=Measure(Ω,10)
    eh=fh-f
    dc=∫(eh*eh)*dΩ
    error_indicators=map(local_views(dc)) do dc
        sqrt.(get_array(dc))
    end
end

function _adapt_model(ranks,model)
    cell_partition=get_cell_gids(model.octree_model.dmodel)
    ref_coarse_flags=map(ranks,partition(cell_partition)) do rank,indices
        flags=zeros(Cint,length(indices))
        flags.=refine_flag
    end
    GridapP4est.adapt(model,ref_coarse_flags)
end

function _adapt_model(ranks,model,error_indicators)
  cell_partition=get_cell_gids(model.octree_model.dmodel)
  ref_coarse_flags=map(ranks,partition(cell_partition)) do rank,indices
      flags=zeros(Cint,length(indices))
      flags.=nothing_flag
  end
  ref_fraction=0.2
  coarsen_fraction=0.05
  adaptive_strategy=FixedFractionAdaptiveFlagsMarkingStrategy(ref_fraction,coarsen_fraction)
  update_adaptivity_flags!(ref_coarse_flags,
                          adaptive_strategy,
                          partition(cell_partition),
                          error_indicators;
                          verbose=true)
  GridapP4est.adapt(model,ref_coarse_flags)
end

function adapt_model(ranks,model,error_indicators;uniform=true)
  if uniform
    println("uniform refinement - refine everywhere")
    _adapt_model(ranks,model)
  else
    println("adaptive refinement - refine on error")
    _adapt_model(ranks,model,error_indicators)
  end
end

function setup_fe_space(model)
    reffe=ReferenceFE(lagrangian,Float64,1)
    Vh=FESpace(model,reffe)
end

models = []
with_mpi() do distribute
    n_refs = 2 # number of refinements
    ranks = distribute(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),)))
    model = CubedSphereDiscreteModel(ranks,n_refs;adaptive=true)
    Vh = setup_fe_space(model)

    fh = interpolate(f,Vh)
    writevtk(Triangulation(model),joinpath(datadir("models"),"cubed_sphere_amr_step_0"),
                cellfields=["f"=>fh],append=false)
    push!(models,model)
    for step=1:n_refs
       fh = interpolate(f,Vh)
       error_indicators=calculate_error_indicators(model,fh)
       model,_ = adapt_model(ranks, model,error_indicators)
        Vh = setup_fe_space(model)
        fh = interpolate(f,Vh)
        writevtk(Triangulation(model),joinpath(datadir("models"),"cubed_sphere_amr_step_$(step)"),
                cellfields=["f"=>fh],append=false)
        push!(models,model)
    end
end

using GridapSolvers
using Gridap.Helpers: @check
import GridapSolvers.MultilevelTools: ModelHierarchyLevel, HierarchicalArray
reverse!(models)


function CubedSphereModelHierarchy(models)
  nlevs = length(models)
  # @check all(map(i -> isa(models[i],GridapDistributed.DistributedAdaptedDiscreteModel),1:nlevs-1)) "Hierarchy models are not adapted models."
  # for lev in 1:nlevs-1
  #   @check Adaptivity.is_child(models[lev],models[lev+1]) "Incorrect hierarchy of models."
  # end
  ranks = get_parts(models[1])
  @check all(m -> length(get_parts(m)) === length(ranks), models) "Models have different communicators."

  level_parts = fill(ranks,nlevs)
  meshes = Vector{ModelHierarchyLevel}(undef,nlevs)
  for lev in 1:nlevs-1
    model = models[lev]
    glue  = model.octree_model.non_conforming_glue #Gridap.Adaptivity.get_adaptivity_glue(models[lev])

    meshes[lev] = ModelHierarchyLevel(lev,model,glue,nothing,nothing)
  end
  meshes[nlevs] = ModelHierarchyLevel(nlevs,models[nlevs],nothing,nothing,nothing)
  return HierarchicalArray(meshes,level_parts)
end


mh = CubedSphereModelHierarchy(models)
