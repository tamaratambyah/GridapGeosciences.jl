"""
In this test, benchmark the assembly of the grad-grad term in the parametric model
vs. the ambient model on a single panel
"""

using GridapGeosciences
using Gridap
using BenchmarkTools
using DrWatson
using Gridap.Geometry
using CountFlops

include("single_panel_ambient_model.jl")
include("analytic_metric.jl")

function bm_func(biform,trial,test)
  assemble_matrix(biform,trial,test)
end



nC = 2^1 ## num cells per panel
radius = 1.0
p_fe = 1
panel_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(nC,nC))
ambient_model = SingleAmbientDiscreteModel(panel_model,1.0)


# function benchmark_gradgrad(ambient_model,panel_model,p_fe,dir,n=nref(ambient_model))
  # println("nref = $n")
  degree = 5*(p_fe+1)

  ################################################################################
  ########## Ambient model
  ################################################################################
  Ω_ambient = Triangulation(ambient_model)
  dΩ_ambient = Measure(Ω_ambient,degree)
  V_ambient = TestFESpace(Ω_ambient, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U_ambient = TrialFESpace(V_ambient)
  grad_grad_ambient(u,v) = ∫( ∇(u)⋅∇(v)   )dΩ_ambient

  @benchmark bm_func($grad_grad_ambient,$U_ambient,$V_ambient)
  t_ambient = @belapsed bm_func($grad_grad_ambient,$U_ambient,$V_ambient)
  flops_ambient = @count_ops bm_func($grad_grad_ambient,$U_ambient,$V_ambient) ### FAILS!


  ################################################################################
  ########## Parametric model -- integrate in the physical domains (includes jacobians)
  ################################################################################
  Ω_panel = Triangulation(panel_model)
  dΩ_panel = Measure(Ω_panel,degree)
  V_panel = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U_panel = TrialFESpace(V_panel)
  grad_grad_panel(u,v) =  ∫( ( ∇(v)⋅ (inv_g⋅ ∇(u) ) )*_sqrtg )dΩ_panel

  @benchmark bm_func($grad_grad_panel,$U_panel,$V_panel)
  t_panel = @belapsed bm_func($grad_grad_panel,$U_panel,$V_panel)
  flops_panel = @count_ops bm_func($grad_grad_panel,$U_panel,$V_panel)


  ################################################################################
  ########## Parametric model -- integrate in reference domain (no jacobians)
  ########## Now the metric goes from the refee -> ambient space
  ################################################################################
  cmap = get_cell_map(panel_model.grid)
  g_ref = lazy_map(Operation(inv_g),cmap)
  g_ref_cf = Gridap.CellData.GenericCellField(g_ref,Ω_panel,ReferenceDomain())

  meas_ref = lazy_map(Operation(_sqrtg),cmap)
  meas_ref_cf = Gridap.CellData.GenericCellField(meas_ref,Ω_panel,ReferenceDomain())

  quad = CellQuadrature(Ω_panel,degree;
          data_domain_style=ReferenceDomain(),
          integration_domain_style=ReferenceDomain())
  dΩ_ref = Gridap.CellData.GenericMeasure(quad)

  # grad_grad_panel_ref(u,v) =  ∫( ( ∇(v)⋅ ∇(u) )*meas_ref_cf )dΩ_ref
  grad_grad_panel_ref(u,v) =  ∫( ( ∇(v)⋅ (g_ref_cf⋅ ∇(u) ) )*meas_ref_cf )dΩ_ref

  @benchmark  bm_func($grad_grad_panel_ref,$U_panel,$V_panel)
  t_ref = @belapsed bm_func($grad_grad_panel_ref,$U_panel,$V_panel)
  flops_ref = @count_ops bm_func($grad_grad_panel_ref,$U_panel,$V_panel)

  ;
  # ### convergence output for DrWatson
  # dir_convergence = dir
  # ( !isdir(dir_convergence)) && mkpath(dir_convergence)
  # output = @strdict n t_panel t_ambient alloc_panel alloc_ambient p_fe  t_panel_ref alloc_panel_ref
  # safesave(datadir(dir_convergence, ("benchmark_nref$(n)_p$(p_fe).jld2")), output)

  # println("Time: Ambient $t_ambient, Parametric $t_panel, Reference $t_panel_ref")
  # println("Allocats: Ambient $alloc_ambient, Parametric $alloc_panel, Reference $alloc_panel_ref")

  # return t_ambient,alloc_ambient, t_panel, alloc_panel



# p_fe = 1
# radius = 1
# n_ref_lvls = 5


# include("single_panel_ambient_model.jl")
# nC = 2^1
#   single_panel_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(nC,nC))
#   single_ambient_model = SingleAmbientDiscreteModel(single_panel_model,radius)

# dir = datadir("grad_grad_single_panel_flops")
# for i in collect(1:n_ref_lvls)
#   nC = 2^i
#   single_panel_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(nC,nC))
#   single_ambient_model = SingleAmbientDiscreteModel(single_panel_model,radius)
#   benchmark_gradgrad(single_ambient_model,single_panel_model,p_fe,dir,i)
# end
