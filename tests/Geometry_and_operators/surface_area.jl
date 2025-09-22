"""
test the integration of surface area
i.e. surface area = ∫ᵧ 1 = ∫ 1 √g
"""

function surface_area(panel_model,degree::Int)
  extact_area = 4*π*RADIUS^2

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)

  meas_cf = CellField(sqrtg,Ω_panel)
  computed_area = sum( ∫( 1.0*meas_cf )dΩ )

  e = abs(computed_area-extact_area)/extact_area
  return e
end

function surface_area_errors(panel_model,degree::Int)
  e = surface_area(panel_model,degree)
  return e,false,false
end


function surface_area_convergence_test(n_ref_lvls)

  for (j,r) in enumerate([1]) #enumerate([1,2,3,4])
    # const global RADIUS = r

    plot()
    for (i,degree) in enumerate([2,4,6,8])
      errs,ns,dxs,slope = convergence_test(surface_area_errors,n_ref_lvls,degree)
      plot_convergence(errs,ns,dxs,slope;leginf=["d=$degree"],colors=[palette(:tab10)[i]])
    end
    savefig(plotsdir()*"/surface_area_convergence_r$j")
  end
end
