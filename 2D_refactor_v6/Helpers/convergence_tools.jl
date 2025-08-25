"""
get_refined_models
returns an array of refined models where
  models[1] == most refined model
  models[end] == coarsest model
"""
function get_refined_models(n_ref_lvls::Int)

  panel_model = coarse_parametric_model()
  panel_models = Vector{DiscreteModel}(undef,n_ref_lvls)

  for n in n_ref_lvls:-1:1
    panel_model = Gridap.Adaptivity.refine(panel_model)
    panel_models[n] = panel_model
  end

  panel_models
end

function h_convergence_test(f,models,fargs...)
  errs = Float64[]
  errs_g = []

  for model in models
    e,eg = f(model,fargs...)
    push!(errs,e)
    push!(errs_g,eg)
  end

  errs, errs_g
end

function convergence_test(f,n_ref_lvls,fargs...)
  models  = get_refined_models(n_ref_lvls)
  errs, errs_g = h_convergence_test(f,models,fargs...)

  ns = map(x->nc(x),models)
  dxs = map(x->dx(nc(x)),models)
  slope = convergence_rate(dxs,errs)

  if typeof(errs_g[1]) == Bool
    return errs,ns,dxs,slope
  else
    return [errs;errs_g],ns,dxs,slope
  end


end

function plot_convergence(errs,ns,dxs,slope;kwargs...)
  r = string(Int(round(slope))) # approximate convergence rate

  plot_error(ns,errs;kwargs...)
  plot!(yscale=:log10,framestyle=:box,
  xscale=:log10,xlabel="n cells",ylabel="L2(u - uh)"
  )
  plot_error(ns,dxs.^slope*errs[1]*1000;leginf=["dx^$r"],colors=kwargs[:colors],ls=[:dash],markers = [:none])
end



## nc = num cells per panel
nc(panel_model) = num_cells(panel_model)/6
dx(nc) = sqrt( 4*π*RADIUS^2 / (6*sqrt(nc)^2) )
nref(nc) = Int(log2(sqrt(nc))) ## level of refinement

function convergence_rate(dxs,errors)
  x = log10.(dxs)
  y = log10.(abs.(errors))
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
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
