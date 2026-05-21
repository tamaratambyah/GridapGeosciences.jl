"""
In this test, benchmark the assembly of the grad-grad term in the parametric model
vs. the ambient model on a single panel
"""

using GridapGeosciences
using Gridap
using BenchmarkTools
using DrWatson
using Gridap.Geometry
# using CountFlops
using GFlops

include("single_panel_ambient_model.jl")
include("analytic_metric.jl")
include("overloads.jl")

function bm_func(biform,trial,test)
  assemble_matrix(biform,trial,test)
end


function total_counts(flops::GFlops.Counter)
  total = 0
  for pn in propertynames(flops)
      total += getproperty(flops, pn)
  end
  return total
end

biform(inv_g,sqrtg,dΩ) = (u,v) -> ∫( ( ∇(v)⋅ (inv_g⋅ ∇(u) ) )*sqrtg )dΩ

function set_up(model,reffee,conf)
  test = TestFESpace(model, reffee; conformity=conf)
  trial = TrialFESpace(test)
  return trial, test
end

nref = 2 ## num cells per panel
radius = 1.0
p_fe = 1
reffee = ReferenceFE(lagrangian,Float64,p_fe)
conf = :H1

function run(nref,dir,reffee,conf,p_fe=1,radius=1.0)
ops_dir = dir*"/operations"
!isdir(ops_dir) && mkpath(ops_dir)

ts_dir = dir*"/time"
!isdir(ts_dir) && mkpath(ts_dir)

flops_dir = dir*"/flops"
!isdir(flops_dir) && mkpath(flops_dir)

nC = 2^nref
panel_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(nC,nC))
ambient_model = SingleAmbientDiscreteModel(panel_model,radius)

degree = 5*(p_fe+1)

################################################################################
########## Ambient model
################################################################################
Ω_ambient = Triangulation(ambient_model)
dΩ_ambient = Measure(Ω_ambient,degree)
U_ambient, V_ambient = set_up(ambient_model,reffee,conf)
grad_grad_ambient(u,v) = ∫( ∇(u)⋅ ∇(v)   )dΩ_ambient
# grad_grad_ambient(u,v) = biform(identity_g,id_sqrtg,dΩ_ambient)(u,v)

bm_func(grad_grad_ambient,U_ambient,V_ambient) # warmup

# @benchmark bm_func($mass_ambient,$U_ambient,$V_ambient)
t_ambient = @belapsed bm_func($grad_grad_ambient,$U_ambient,$V_ambient)
ops_ambient = @count_ops bm_func($grad_grad_ambient,$U_ambient,$V_ambient)
flops_ambient = @gflops bm_func($grad_grad_ambient,$U_ambient,$V_ambient)



################################################################################
########## Parametric model -- integrate in the physical domain (includes jacobians)
################################################################################
Ω_panel = Triangulation(panel_model)
dΩ_panel = Measure(Ω_panel,degree)
U_panel, V_panel = set_up(panel_model,reffee,conf)
# grad_grad_panel(u,v) =  ∫( ( ∇(v)⋅ (inv_g⋅ ∇(u) ) )*_sqrtg )dΩ_panel
grad_grad_panel(u,v) =  biform(inv_g,_sqrtg,dΩ_panel)(u,v)

bm_func(grad_grad_panel,U_panel,V_panel) # warmup

# @benchmark bm_func($mass_panel,$U_panel,$V_panel)
t_panel = @belapsed bm_func($grad_grad_panel,$U_panel,$V_panel)
ops_panel = @count_ops bm_func($grad_grad_panel,$U_panel,$V_panel)
flops_panel = @gflops bm_func($grad_grad_panel,$U_panel,$V_panel)



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

# grad_grad_panel_ref(u,v) =  ∫( ( ∇(v)⋅ (g_ref_cf⋅ ∇(u) ) )*meas_ref_cf )dΩ_ref
grad_grad_panel_ref(u,v) =  biform(g_ref_cf,meas_ref_cf,dΩ_ref)(u,v)

bm_func(grad_grad_panel_ref,U_panel,V_panel) # warmup

# @benchmark  bm_func($mass_panel_ref,$U_panel,$V_panel)
t_ref = @belapsed bm_func($grad_grad_panel_ref,$U_panel,$V_panel)
ops_ref = @count_ops bm_func($grad_grad_panel_ref,$U_panel,$V_panel)
flops_ref = @gflops bm_func($grad_grad_panel_ref,$U_panel,$V_panel)

total_ops_ambient = total_counts(ops_ambient)
total_ops_panel = total_counts(ops_panel)
total_ops_ref = total_counts(ops_ref)


keys = ["ambient","panel","ref","nref"]

ops_vals = [total_ops_ambient, total_ops_panel, total_ops_ref, nref]
ops_dict = Dict(keys .=> ops_vals)
safesave(datadir(ops_dir, ("ops_n$(nref)_p$p_fe.jld2")), ops_dict)

ts_vals = [t_ambient, t_panel, t_ref, nref]
ts_dict = Dict(keys .=> ts_vals)
safesave(datadir(ts_dir, ("ts_n$(nref)_p$p_fe.jld2")), ts_dict)

flops_vals = [flops_ambient, flops_panel, flops_ref, nref]
flops_dict = Dict(keys .=> flops_vals)
safesave(datadir(flops_dir, ("flops_n$(nref)_p$p_fe.jld2")), flops_dict)


end


dir = datadir("gradgrad")
!isdir(dir) && mkpath(dir)
for nref in collect(1:5)
  run(nref,dir,reffee,conf)
end
