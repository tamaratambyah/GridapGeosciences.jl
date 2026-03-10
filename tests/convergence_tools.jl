using Gridap
using Gridap.Adaptivity, Gridap.Helpers
using GridapDistributed
using GridapGeosciences
using PartitionedArrays
using Plots, LaTeXStrings
using JLD2
using Test

include("helpers.jl")

function foldername(name,octree=false,threedims=false)
  dir = datadir(name)

  if threedims
    return dir*"_3D"
  end

  if octree
    return dir*"_Octree"
  end

  return dir
end


function get_models(ranks,nprocs,n_ref_lvls::Int;threedims=false,octree=false)
  s_models  = get_refined_models(n_ref_lvls)

  if threedims
    i_am_main(ranks) && println("3D models")
    return get_3D_octree_refined_models(ranks,n_ref_lvls)
  end

  if octree
    i_am_main(ranks) && println("Octree models")
    return get_octree_refined_models(ranks,n_ref_lvls)
  end

  if nprocs > 1
    i_am_main(ranks) && println("Distributed test")
    models,  = get_distributed_refined_models(ranks,nprocs,s_models)
    return models
  else
    i_am_main(ranks) && println("Serial test")
    return s_models
  end

end


"""
get_refined_models
returns an array of refined models where
  models[1] == most refined model
  models[end] == coarsest model
"""
function get_refined_models(n_ref_lvls::Int,coarse_model=false)

  panel_model = coarse_parametric_model()
  panel_models = Vector{DiscreteModel}(undef,n_ref_lvls)

  for n in n_ref_lvls:-1:1
    panel_model = Gridap.Adaptivity.refine(panel_model)
    panel_models[n] = panel_model
  end

  if coarse_model
    push!(panel_models,coarse_parametric_model())
  end

  panel_models
end

function get_distributed_refined_models(ranks,nprocs,n_ref_lvls::Int,coarse_s_model=false)
  s_models  = get_refined_models(n_ref_lvls,coarse_s_model)
  dmodels, dpanel_ids, owned_panel_ids = get_distributed_refined_models(ranks,nprocs,s_models)
  dmodels
end


function get_distributed_refined_models(ranks,nprocs,s_models::Vector{<:DiscreteModel})

  # get refined models in serial
  spanel_ids = map(m->get_panel_ids(m),s_models)
  s_model_coarse = s_models[end]

  # extract the models and glues in arrays
  models = map(m->Adaptivity.get_model(m),s_models[1:end-1])
  glues = map(m->get_adaptivity_glue(m),s_models[1:end-1])
  if typeof(s_model_coarse) <: AdaptedDiscreteModel
    push!(models,Adaptivity.get_model(s_model_coarse))
  else
    push!(models,s_model_coarse)
  end

  # partition the processors
  part_to_cells = [PartitionedArrays.local_range(rank,nprocs,num_cells(s_model_coarse)) for rank in 1:nprocs]

  # store the partition for each model
  cell_to_part = Vector{Vector{Int32}}(undef,length(models))

  # get the coarse partition
  coarse_cell_to_part = zeros(Int32,num_cells(s_model_coarse))
  for (rank, cells) in enumerate(part_to_cells)
    coarse_cell_to_part[cells] .= rank
  end
  cell_to_part[end] = coarse_cell_to_part

  # get the refine partition based on glue
  for level in length(models)-1:-1:1
    n2o_cells = glues[level].n2o_faces_map[3]
    cell_to_part[level] = cell_to_part[level+1][n2o_cells]
  end

  # construct array of distributed models, distributed panel ids (all + owned only)
  dmodels = Vector{DistributedParametricDiscreteModel}(undef,length(models))
  dpanel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))
  owned_panel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))

  # the coarsest model is the last in the list
  coarse_dmodel = DiscreteModel(ranks,models[end],cell_to_part[end])
  dpanel_ids[end],owned_panel_ids[end] = distributed_panel_ids(coarse_dmodel,spanel_ids[end])
  dmodels[end] = DistributedParametricDiscreteModel(coarse_dmodel,dpanel_ids[end])

  # loop backwards through refinement levels
  for level in length(models)-1:-1:1
    child = DiscreteModel(ranks,models[level],cell_to_part[level])
    parent = dmodels[level+1]
    glue = DistributedAdaptivityGlue(glues[level],parent,child)
    dpanel_ids[level],owned_panel_ids[level] = distributed_panel_ids(child,spanel_ids[level])

    dmodels[level] = DistributedParametricDiscreteModel(child, dpanel_ids[level])
  end

  return dmodels, dpanel_ids, owned_panel_ids

end

#### get array of dmodel - octree models
function get_octree_refined_models(ranks,n_ref_lvls::Int,coarse_model=false)

  dmodels = Vector{DistributedParametricDiscreteModel}(undef,n_ref_lvls)

  for (i,n) in enumerate(n_ref_lvls:-1:1)
    parametric_octree_dmodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=n)
    dmodels[i] = parametric_octree_dmodel.parametric_dmodel
  end


  if coarse_model
    parametric_octree_dmodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=0)
    push!(dmodels,parametric_octree_dmodel.parametric_dmodel)
  end

  return dmodels
end

#### get array of octree models
function _get_octree_refined_models(ranks,n_ref_lvls::Int,coarse_model=false)

  omodels = Vector{ParametricOctreeDistributedDiscreteModel}(undef,n_ref_lvls)

  for (i,n) in enumerate(n_ref_lvls:-1:1)
    parametric_octree_model = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=n)
    omodels[i] = parametric_octree_model
  end


  if coarse_model
    parametric_octree_model = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=0)
    push!(omodels,parametric_octree_model)
  end

  return omodels
end

function get_3D_octree_refined_models(ranks,n_ref_lvls::Int)

  dmodels = Vector{DistributedParametricDiscreteModel}(undef,n_ref_lvls)

  for (i,n) in enumerate(n_ref_lvls:-1:1)
    octree3_model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
                        num_horizontal_uniform_refinements=n,
                        num_vertical_uniform_refinements=n);
    dmodels[i] = octree3_model.parametric_dmodel
  end

  return dmodels

end

function get_3D_octree_horizontal_refined_models(ranks,n_ref_lvls_horiztontal::Int,num_vertical_uniform_refinements=3)
  n_ref_lvls = n_ref_lvls_horiztontal

  dmodels = Vector{DistributedParametricDiscreteModel}(undef,n_ref_lvls)

  for (i,n) in enumerate(n_ref_lvls:-1:2)
    octree3_model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
                        num_horizontal_uniform_refinements=n,
                        num_vertical_uniform_refinements=num_vertical_uniform_refinements);
    dmodels[i] = octree3_model.parametric_dmodel
  end

  return dmodels

end

function get_3D_octree_vertical_refined_models(ranks,num_vertical_uniform_refinements::Int,n_ref_lvls_horiztontal=3)
  n_ref_lvls = num_vertical_uniform_refinements

  dmodels = Vector{DistributedParametricDiscreteModel}(undef,n_ref_lvls)

  for (i,n) in enumerate(n_ref_lvls:-1:1)
    octree3_model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
                        num_horizontal_uniform_refinements=n_ref_lvls_horiztontal,
                        num_vertical_uniform_refinements=n);
    dmodels[i] = octree3_model.parametric_dmodel
  end

  return dmodels

end


function get_ranks(model::DiscreteModel)
    return [true]
end

function get_ranks(dmodel::GridapDistributed.DistributedDiscreteModel)
  get_parts(dmodel)
end


"""
convergence test functions
"""

function _h_convergence_test(models,f,p_fe::Int,dir::String,fargs...)
  errs = Float64[]
  errs_g = []
  errs_f = []

  for model in models
    e,eg,ef = f(model,p_fe,dir,fargs...)
    push!(errs,e)
    push!(errs_g,eg)
    push!(errs_f,ef)
  end

  errs, errs_g, errs_f
end

function h_convergence_test(models::AbstractArray,f,p_fe::Int,dir::String,fargs...)
  errs, errs_g, errs_f = _h_convergence_test(models,f,p_fe,dir,fargs...)

  ns = map(x->nc(x),models)
  dxs = map(x->dx(x),models)
  slope = convergence_rate(dxs,errs)

  if typeof(errs_g[1]) == Bool
    return errs,ns,dxs,slope
  elseif typeof(errs_f[1]) == Bool
    return [errs;errs_g],ns,dxs,slope
  else
    return [errs;errs_g;errs_f],ns,dxs,slope
  end
end

function h_convergence_test_vertical(models::AbstractArray,f,p_fe::Int,dir::String,fargs...)
  errs, errs_g, errs_f = _h_convergence_test(models,f,p_fe,dir,fargs...)

  ns = map(x->_nc_vertical(x),models)
  dxs = map(x->1/_nc_vertical(x),models)
  slope = convergence_rate(dxs,errs)

  if typeof(errs_g[1]) == Bool
    return errs,ns,dxs,slope
  elseif typeof(errs_f[1]) == Bool
    return [errs;errs_g],ns,dxs,slope
  else
    return [errs;errs_g;errs_f],ns,dxs,slope
  end
end

# set ranks = [true] for serial
function p_convergence_test(ranks,ps::Vector{Int},models::AbstractArray,convergence_func,dir::String,fargs...)
  # i_am_main(ranks) && println("auto convergence test")

  errors = Vector{Vector{Float64}}(undef,length(ps))
  ns = Vector{Vector{Float64}}(undef,length(ps))
  dxs = Vector{Vector{Float64}}(undef,length(ps))
  slopes = Vector{Float64}(undef,length(ps))


  for (i,p_fe) in enumerate(ps)
    # i_am_main(ranks) && println("p_fe = $p_fe")
    errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,convergence_func,p_fe,dir,fargs...)
    i_am_main(ranks) && print_convergence_results(errors[i],ns[i],dxs[i],slopes[i],p_fe)

    output = @strdict errors ns dxs slopes
    i_am_main(ranks) && safesave(datadir(dir, ("convergence_p$p_fe.jld2")), output)

    # @test slopes[i] >= p_fe + 1
  end

  output = @strdict errors ns dxs slopes ps
  i_am_main(ranks) && safesave(datadir(dir, ("convergence.jld2")), output)

end

function convergence_rate(dxs,errors)
  x = log10.(dxs)
  y = log10.(abs.(errors))
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end


"""
plotting and printing helper functions
"""

function plot_convergence(errs,ns,dxs,slope;kwargs...)
  r = string(Int(round(slope))) # approximate convergence rate

  plot_error(ns,errs;kwargs...)
  plot!(yscale=:log10,framestyle=:box,
  xscale=:log10,xlabel="n cells",ylabel="L2(u - uh)"
  )
  plot_error(ns,dxs.^slope*errs[1]*1000;leginf=["dx^$r"],colors=kwargs[:colors],ls=[:dash],markers = [:none])
end


function plot_error(ns,errs;
  leginf = fill(false,Int(length(errs)/length(ns))),
  ls=[:solid, :dash, :dot, :dashdot, :dashdotdot],
  colors = palette(:tab10),
  markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross],
  ms=[6,6,8,6,6,8,8] )

  nsims = Int(length(errs)/length(ns))
  for i in 1:nsims
    idx1 = 1 + (i-1)*length(ns)
    idx2 = (i)*length(ns)
    plot!(ns,
          errs[idx1:idx2],
          lw=3,
          markersize=ms[i],
          c=colors[i],ls=ls[i], markershape=markers[i],
          label=leginf[i])
  end

end


function plot_convergence_from_saved(dir,simName,varNames=["u"])
  dd = load(datadir(dir, ("$simName.jld2")))
  dxs = dd["dxs"]
  errors = dd["errors"]
  ns = dd["ns"]
  ps = dd["ps"]
  slopes = dd["slopes"]


  cidx = ps
  if any(iszero.(ps))
    cidx = ps .+ 1
  end

  plot()
  for (i,p_fe) in enumerate(ps)
    leginf = map(x->"$x: p=$p_fe", varNames)
    cols = map(x->palette(:tab10)[cidx[i]], varNames )
    plot_convergence(errors[i],ns[i],dxs[i],slopes[i];leginf=leginf,colors=cols)
  end
  savefig(dir*"/$simName")
end


function print_convergence_results(errors::AbstractArray,ns::AbstractArray{Vector{Float64}},
  dxs::AbstractArray{Vector{Float64}},slopes::AbstractArray{Float64},ps::AbstractArray{Int})
  for (i,p_fe) in enumerate(ps)
    print_convergence_results(errors[i],ns[i],dxs[i],slopes[i],p_fe)
  end
end

function print_convergence_results(errors::AbstractArray,ns::AbstractArray{Float64},dxs::AbstractArray{Float64},slopes::Float64,p_fe::Int)
    println("p_fe = $p_fe, slope = ", slopes)
    println("\t Error: ", errors)
    println("\t    nc: ", ns)
    println("\t    dx: ", dxs)
end
