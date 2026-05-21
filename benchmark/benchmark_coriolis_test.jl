"""
In this test, benchmark the assembly of the coriolis term in the parametric model
vs. the ambient model
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

Aperp = [0 -1
        1 0]
Rperp = TensorValue(Aperp)
biform(dΩ) = (u,v) ->  ∫( ( ( (Rperp⋅ u)⋅v))  )dΩ

function set_up(model,reffee,conf)
  test = TestFESpace(model, reffee; conformity=conf)
  trial = TrialFESpace(test)
  return trial, test
end

nref = 2 ## num cells per panel
radius = 1.0
p_fe = 1
reffee = ReferenceFE(raviart_thomas,Float64,p_fe)
conf = :HDiv

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

coriolis_term_ambient(u,v) = ∫( ( ( normal_vec × u)⋅v)  )dΩ_ambient

bm_func(coriolis_term_ambient,U_ambient,V_ambient) # warmup

# @benchmark bm_func($mass_ambient,$U_ambient,$V_ambient)
t_ambient = @belapsed bm_func($coriolis_term_ambient,$U_ambient,$V_ambient)
ops_ambient = @count_ops bm_func($coriolis_term_ambient,$U_ambient,$V_ambient)
flops_ambient = @gflops bm_func($coriolis_term_ambient,$U_ambient,$V_ambient)



################################################################################
########## Parametric model -- integrate in the physical domain (includes jacobians)
################################################################################
Ω_panel = Triangulation(panel_model)
dΩ_panel = Measure(Ω_panel,degree)
U_panel, V_panel = set_up(panel_model,reffee,conf)


coriolis_term_panel(u,v) =  biform(dΩ_panel)(u,v)

bm_func(coriolis_term_panel,U_panel,V_panel) # warmup

# @benchmark bm_func($coriolis_term_panel,$U_panel,$V_panel)
t_panel = @belapsed bm_func($coriolis_term_panel,$U_panel,$V_panel)
ops_panel = @count_ops bm_func($coriolis_term_panel,$U_panel,$V_panel)
flops_panel = @gflops bm_func($coriolis_term_panel,$U_panel,$V_panel)



################################################################################
########## Parametric model -- integrate in reference domain (no jacobians)
########## Now the metric goes from the refee -> ambient space
################################################################################
quad = CellQuadrature(Ω_panel,degree;
        data_domain_style=ReferenceDomain(),
        integration_domain_style=ReferenceDomain())
dΩ_ref = Gridap.CellData.GenericMeasure(quad)

coriolis_term_panel_ref(u,v) =  biform(dΩ_ref)(u,v)

bm_func(coriolis_term_panel_ref,U_panel,V_panel) # warmup

# @benchmark  bm_func($mass_panel_ref,$U_panel,$V_panel)
t_ref = @belapsed bm_func($coriolis_term_panel_ref,$U_panel,$V_panel)
ops_ref = @count_ops bm_func($coriolis_term_panel_ref,$U_panel,$V_panel)
flops_ref = @gflops bm_func($coriolis_term_panel_ref,$U_panel,$V_panel)

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


dir = datadir("coriolis")
!isdir(dir) && mkpath(dir)
for nref in collect(1:5)
  run(nref,dir,reffee,conf)
end
