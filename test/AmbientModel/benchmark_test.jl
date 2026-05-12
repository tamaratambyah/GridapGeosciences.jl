"""
In this test, benchmark the assembly of the coriolis term in the parametric model
vs. the ambient model

Output: bm_ambient: 11.1 ms, mem 6.38 MiB, alloc 79759
        bm_panel:   2.66 ms, mem 2.70 MiB, alloc 20554
Conclusion: Parametric model is much better for performance for coriolis
most likely due to cross product.
"""

using GridapGeosciences
using Gridap
using BenchmarkTools


function get_setup(model,p_fe)
  degree = 5*(p_fe+1)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  V = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  return Ω,dΩ,U,V
end

function fcor_X(xyz)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    sin(ϕ)
end


function benchmark_coriolis(ambient_model,p_fe)

  ################################################################################
  ########## Ambient model
  ################################################################################
  Ω_ambient,dΩ_ambient,U_ambient,V_ambient = get_setup(ambient_model,p_fe)
  cor_cf_ambient = CellField(fcor_X,Ω_ambient)

  ## Here we construct the coriolis term on the surface: ∫( ̃f ( ̃k × ̃u  )  )dΩ
  n_surf = get_surface_normal(Ω_ambient)
  coriolis_term_ambient(u,v) = ∫( cor_cf_ambient*( ( n_surf × u)⋅v)  )dΩ_ambient

  bm_ambient() = assemble_matrix(coriolis_term_ambient,U_ambient,V_ambient)



  ################################################################################
  ########## Parametric model
  ################################################################################
  panel_model = get_parametric_model(ambient_model)

  Ω_panel,dΩ_panel,U_panel,V_panel = get_setup(panel_model,p_fe)

  cor_cf_panel = ParametricCellField(panel_to_cartesian(fcor_X),Ω_panel)

  ## Here we construct the coriolis term in 2D using the rotation matrx
  Aperp = [0 -1
          1 0]
  Rperp = TensorValue(Aperp)
  Rperp_cf = CellField(Rperp,Ω_panel)

  coriolis_term_panel(u,v) = ∫( ( cor_cf_panel*( (Rperp_cf⋅ u)⋅v))  )dΩ_panel

  bm_panel() = assemble_matrix(coriolis_term_panel,U_panel,V_panel)


  # @benchmark bm_ambient()
  # @benchmark bm_panel()

  ## elapsed time
  t_ambient = @belapsed bm_ambient()
  t_panel = @belapsed bm_panel()

  ## allocations
  alloc_ambient = @ballocations bm_ambient()
  alloc_panel = @ballocations bm_panel()

  return t_ambient,alloc_ambient, t_panel, alloc_panel

end

p_fe = 1
radius = 1
n_ref_lvls = 4
ambient_models = get_ambient_refined_models(n_ref_lvls,radius)

ts_am = Vector{Float64}(undef,n_ref_lvls)
ts_pm = Vector{Float64}(undef,n_ref_lvls)

alcs_am = Vector{Float64}(undef,n_ref_lvls)
alcs_pm = Vector{Float64}(undef,n_ref_lvls)

for (i,model) in enumerate(ambient_models)
  ts_am[i],alcs_am[i], ts_pm[i], alcs_pm[i] = benchmark_coriolis(model,p_fe)
end

t_am2pm = ts_am./ts_pm
a_am2pm = alcs_am./alcs_pm
lvs = map(x->nref(x),ambient_models)

dir = @__DIR__
using Plots
plot()
plot!(lvs,t_am2pm,lw=2,label="time")
plot!(lvs,a_am2pm,lw=2,label="allocations")
plot!(xlabel="lvl",ylabel="ambient/panel")
xl = map(x->string(Int(((x)))),lvs)
plot!(xticks = (lvs, xl))
savefig(dir*"/benchmark.png")
