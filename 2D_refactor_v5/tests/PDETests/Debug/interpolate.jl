using DrWatson
using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots, LaTeXStrings
include("../../../src/initialise.jl")
include("../pde_helpers.jl")


function uX_scalar(X)
  X[3]^2
end

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
models = get_refined_models(manifold_model,3)

p_fe = 2
degree = 4*(p_fe+1)

errs = []
errs_g = []
errs_amb = []

for manifold_model in models[1:2]
  println("N = ", num_cells(manifold_model))
  panel_ids = get_panel_ids(manifold_model)

  ################################################################################
  ##### ambient space
  ################################################################################
  ambient_model = get_ambient_model(manifold_model)
  Ω_amb = Triangulation(ambient_model)
  dΩ_amb = Measure(Ω_amb,degree)

  cell_field = map(p->GenericField(uX_scalar),panel_ids)
  ucf_ambient = CellData.GenericCellField(cell_field,Ω_amb,PhysicalDomain())

  ################################################################################
  ##### parametric space
  ################################################################################
  ####### check compatibility in parametric space
  Ω = Triangulation(manifold_model)
  m = Metric(cubedsphere,Ω)

  dΩ = Measure(Ω,degree)
  dΩg = Measure(m,Ω,degree)

  cell_field = map(p->GenericField(u_scalar_ambient2parametric(p,uX_scalar)),panel_ids)
  ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

  rhs = ucf
  compat = sum( ∫(rhs)dΩ   )
  println("Compatibility: ", compat)

  V = TestFESpace(manifold_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U = TrialFESpace(V)

  uh = interpolate(rhs,U)

  ## map parametric FEFunction back to ambient space
  cf_ambient = parametric_cf_2_ambient(manifold_model,uh.cell_field)


  #### Compute errors
  e = l2(uh-ucf,dΩ)
  eg = l2(uh-ucf,dΩg)
  e_amb = l2(cf_ambient-ucf_ambient,dΩ_amb)

  push!(errs,e)
  push!(errs_g,eg)
  push!(errs_amb,e_amb)

end

nc = [1,2,4,8]
dx =   ( sqrt.( 4*π*RADIUS^2 ./ (6*nc.^2) ) )

plot()
plot_error(nc,errs)
plot_error(nc,errs_g;colors=[:orange])
plot_error(nc,errs_amb;colors=[:green],ls=[:dashdot])
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:bottomright)
plot!(nc,1e-4dx.^(5.5),lw=2,c=:black,label="dx^5.5")
savefig(plotsdir()*"/interpolate_cubedsphere")
# savefig(plotsdir()*"/L2projection_cubedsphere_williamson")
