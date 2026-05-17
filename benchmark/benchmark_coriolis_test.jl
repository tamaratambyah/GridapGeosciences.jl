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
using DrWatson
using DataFrames


function get_setup(model,p_fe)
  degree = 5*(p_fe+1)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  V = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  return Ω,dΩ,U,V
end


function benchmark_coriolis(ambient_model,panel_model,p_fe,dir)

  ################################################################################
  ########## Ambient model
  ################################################################################
  Ω_ambient,dΩ_ambient,U_ambient,V_ambient = get_setup(ambient_model,p_fe)

  ## Here we construct the coriolis term on the surface: ∫( ̃f ( ̃k × ̃u  )  )dΩ
  n_surf = get_surface_normal(Ω_ambient)
  coriolis_term_ambient(u,v) = ∫( ( ( n_surf × u)⋅v)  )dΩ_ambient

  t_ambient = @belapsed assemble_matrix($coriolis_term_ambient,$U_ambient,$V_ambient)
  alloc_ambient = @ballocations assemble_matrix($coriolis_term_ambient,$U_ambient,$V_ambient)


  ################################################################################
  ########## Parametric model
  ################################################################################
  Ω_panel,dΩ_panel,U_panel,V_panel = get_setup(panel_model,p_fe)

  ## Here we construct the coriolis term in 2D using the rotation matrx
  Aperp = [0 -1
          1 0]
  Rperp = TensorValue(Aperp)
  Rperp_cf = CellField(Rperp,Ω_panel)

  coriolis_term_panel(u,v) = ∫( ( ( (Rperp_cf⋅ u)⋅v))  )dΩ_panel

  bm_panel() = assemble_matrix(coriolis_term_panel,U_panel,V_panel)

  t_panel = @belapsed assemble_matrix($coriolis_term_panel,$U_panel,$V_panel)
  alloc_panel = @ballocations assemble_matrix($coriolis_term_panel,$U_panel,$V_panel)


  ### convergence output for DrWatson
  dir_convergence = dir
  ( !isdir(dir_convergence)) && mkpath(dir_convergence)
  n = nref(ambient_model)
  output = @strdict n t_panel t_ambient alloc_panel alloc_ambient p_fe
  safesave(datadir(dir_convergence, ("benchmark_nref$(n)_p$(p_fe).jld2")), output)


  return t_ambient,alloc_ambient, t_panel, alloc_panel

end

p_fe = 1
radius = 1
n_ref_lvls = 5
ambient_models = get_ambient_refined_models(n_ref_lvls,radius)

ts_am = Vector{Float64}(undef,n_ref_lvls)
ts_pm = Vector{Float64}(undef,n_ref_lvls)

alcs_am = Vector{Float64}(undef,n_ref_lvls)
alcs_pm = Vector{Float64}(undef,n_ref_lvls)

dir = datadir("coriolis")
for (i,model) in enumerate(ambient_models)
  panel_model = get_parametric_model(model)
  ts_am[i],alcs_am[i], ts_pm[i], alcs_pm[i] = benchmark_coriolis(model,panel_model,p_fe,dir)
end
