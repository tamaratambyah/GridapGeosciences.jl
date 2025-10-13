include("Williamson_functions_v2.jl")

function williamson2_convergence_test(ranks::AbstractArray,nprocs::Int,
  solver_errors,ζs=[0.0],n_ref_lvls=4,ps=[1],ls=LUSolver(),return_vtk=false)

  solverName = string(solver_errors)[1:end-7]

  # serial models
  models  = get_refined_models(n_ref_lvls)

  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    models,  = get_distributed_refined_models(ranks,nprocs,models)
  end

  dir = datadir("Williamson2ConvergenceTest_$solverName")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  for (i,ζ) in enumerate(ζs)
    simName =  "williamson2_$(solverName)_convergence_func_z$i"
    i_am_main(ranks) && println(simName)

    h = panel_to_cartesian(h₀(ζ))
    vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
    f = panel_to_cartesian(f₀(ζ))
    η = panel_to_cartesian(η₀(ζ))

    errors = Vector{Vector{Float64}}(undef,length(ps))
    ns = Vector{Vector{Float64}}(undef,length(ps))
    dxs = Vector{Vector{Float64}}(undef,length(ps))
    slopes = Vector{Float64}(undef,length(ps))

    for (i,p_fe) in enumerate(ps)
      i_am_main(ranks) && println("p_fe = $p_fe")
      errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,solver_errors,p_fe,dir,h,vX,f,η,ls,return_vtk)
    end

    i_am_main(ranks) && print_convergence_results(errors,ns,dxs,slopes,ps)

    output = @strdict errors ns dxs slopes ps
    i_am_main(ranks) && safesave(datadir(dir, ("$simName.jld2")), output)

    i_am_main(ranks) && plot_convergence_from_saved(dir,simName,["u","ϕ"])

  end


end
