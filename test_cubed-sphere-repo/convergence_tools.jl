using Gridap
using Gridap.Adaptivity, Gridap.Helpers
using GridapDistributed
using GridapGeosciences
using PartitionedArrays
using Plots, LaTeXStrings
using JLD2
using Test
using DrWatson

# include("helpers.jl")

"""
get_ranks
"""
function get_ranks(model::DiscreteModel)
  return [true]
end

function get_ranks(dmodel::GridapDistributed.DistributedDiscreteModel)
get_parts(dmodel)
end



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
    return get_octree_refined_models(ranks,n_ref_lvls,radius)
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






#### get array of octree models
function _get_octree_refined_models(ranks,n_ref_lvls::Int,radius,coarse_model=false)

  omodels = Vector{CubedSphere2DParametricOctreeDistributedDiscreteModel}(undef,n_ref_lvls)

  for (i,n) in enumerate(n_ref_lvls:-1:1)
    parametric_octree_model = CubedSphere2DParametricOctreeDistributedDiscreteModel(ranks, radius; num_initial_uniform_refinements=n)
    omodels[i] = parametric_octree_model
  end


  if coarse_model
    parametric_octree_model = CubedSphere2DParametricOctreeDistributedDiscreteModel(ranks, radius; num_initial_uniform_refinements=0)
    push!(omodels,parametric_octree_model)
  end

  return omodels
end



function get_3D_octree_horizontal_refined_models(ranks,n_ref_lvls_horiztontal::Int,num_vertical_uniform_refinements=3)
  n_ref_lvls = n_ref_lvls_horiztontal

  dmodels = Vector{CubedSphere2DParametricDistributedDiscreteModel}(undef,n_ref_lvls)

  for (i,n) in enumerate(n_ref_lvls:-1:2)
    octree3_model = GridapGeosciences.Distributed.CubedSphere3DParametricOctreeDistributedDiscreteModel(
                        ranks,radius,thickness;
                        num_horizontal_uniform_refinements=n,
                        num_vertical_uniform_refinements=num_vertical_uniform_refinements);
    dmodels[i] = octree3_model.parametric_dmodel
  end

  return dmodels

end

function get_3D_octree_vertical_refined_models(ranks,num_vertical_uniform_refinements::Int,n_ref_lvls_horiztontal=3)
  n_ref_lvls = num_vertical_uniform_refinements

  dmodels = Vector{CubedSphere2DParametricDistributedDiscreteModel}(undef,n_ref_lvls)

  for (i,n) in enumerate(n_ref_lvls:-1:1)
    octree3_model = CubedSphere3DParametricOctreeDistributedDiscreteModel(ranks,radius,thickness;
                        num_horizontal_uniform_refinements=n_ref_lvls_horiztontal,
                        num_vertical_uniform_refinements=n);
    dmodels[i] = octree3_model.parametric_dmodel
  end

  return dmodels

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
