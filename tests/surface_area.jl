## make models
cube_model = coarse_cube_model(π/4,6)
cube_model = Gridap.Adaptivity.refine(cube_model)
cube_model = Gridap.Adaptivity.refine(cube_model)

panel_model = parametric_model(cube_model)

function surface_area(panel_model::ParametricDiscreteModel,degree::Int)
  extact_area = 4*π*RADIUS^2

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)

  meas_cf = CellField(sqrtg,Ω_panel)
  computed_area = sum( ∫( 1.0*meas_cf )dΩ )

  e = abs(computed_area-extact_area)/extact_area
  return e,false,false
end

function surface_area_convergence_test(n_ref_lvls)

  for (j,r) in enumerate([1,2,3,4])
    global RADIUS = r

    plot()
    for (i,degree) in enumerate([2,4,6,8])
      errs,ns,dxs,slope = convergence_test(surface_area,n_ref_lvls,degree)
      plot_convergence(errs,ns,dxs,slope;leginf=["d=$degree"],colors=[palette(:tab10)[i]])
    end
    savefig(plotsdir()*"/surface_area_convergence_r$j")
  end
end

surface_area_convergence_test(4)
