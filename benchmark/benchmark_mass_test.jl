"""
In this test, benchmark the assembly of the scalar mass matrix in the parametric model
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

function bm_func(biform,trial,test)
  assemble_matrix(biform,trial,test)
end

## We have to overwrite this function for CountFlops to be compatible with
## UnstructuredDiscreteModel. This is because it does not know that evaluate on
## the iterator is type stable (and for some reason does a ccal).
## In this function, provide the type to vals explicitily.
## The length of this compressed array is one, so it should skip the for loop
# function Gridap.Arrays._lazy_map_compressed(g::CompressedArray...)
#   #vals = map(evaluate, map(gi->gi.values,g)...)
#   n = length(first(g).values)
#   v1 = evaluate(map(gi -> first(gi.values),g)...)
#   vals = Vector{typeof(v1)}(undef,n)
#   vals[1] = v1
#   for i in 2:n
#     vals[i] = evaluate(map(gi -> gi.values[i],g)...)
#   end
#   ptrs = first(g).ptrs
#   Gridap.Arrays.CompressedArray(vals,ptrs)
# end
function Gridap.Arrays._lazy_map_compressed(g::CompressedArray...)
  #vals = map(evaluate, map(gi->gi.values,g)...)
  n = length(first(g).values)
  v1 = evaluate(map(gi -> first(gi.values),g)...)
  vals = Vector{typeof(v1)}(undef,n)
  for i in 1:n
    vals[i] = evaluate(map(gi -> gi.values[i],g)...)
  end
  ptrs = first(g).ptrs
  Gridap.Arrays.CompressedArray(vals,ptrs)
end


function total_counts(flops::GFlops.Counter)
  total = 0
  for pn in propertynames(flops)
      total += getproperty(flops, pn)
  end
  return total
end

biform(sqrtg,dΩ) = (u,v) -> ∫( ( u*v )*sqrtg )dΩ


nC = 2^2 ## num cells per panel
radius = 1.0
p_fe = 1
_panel_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(nC,nC))

panel_model = _panel_model
ambient_model = SingleAmbientDiscreteModel(panel_model,1.0)


degree = 5*(p_fe+1)

################################################################################
########## Ambient model
################################################################################
Ω_ambient = Triangulation(ambient_model)
dΩ_ambient = Measure(Ω_ambient,degree)
V_ambient = TestFESpace(Ω_ambient, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
U_ambient = TrialFESpace(V_ambient)
# mass_ambient(u,v) = ∫( u*v   )dΩ_ambient
mass_ambient(u,v) = biform(id_sqrtg,dΩ_ambient)(u,v)


bm_func(mass_ambient,U_ambient,V_ambient) # warmup

# @benchmark bm_func($mass_ambient,$U_ambient,$V_ambient)
# t_ambient = @belapsed bm_func($mass_ambient,$U_ambient,$V_ambient)
ops_ambient = @count_ops bm_func($mass_ambient,$U_ambient,$V_ambient)
flops_ambient = @gflops bm_func($mass_ambient,$U_ambient,$V_ambient)


################################################################################
########## Parametric model -- integrate in the physical domain (includes jacobians)
################################################################################
Ω_panel = Triangulation(panel_model)
dΩ_panel = Measure(Ω_panel,degree)
V_panel = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
U_panel = TrialFESpace(V_panel)
# mass_panel(u,v) =  ∫( ( u*v )*_sqrtg )dΩ_panel
mass_panel(u,v) =  biform(_sqrtg,dΩ_panel)(u,v)

bm_func(mass_panel,U_panel,V_panel) # warmup

# @benchmark bm_func($mass_panel,$U_panel,$V_panel)
# t_panel = @belapsed bm_func($mass_panel,$U_panel,$V_panel)
ops_panel = @count_ops bm_func($mass_panel,$U_panel,$V_panel)
flops_panel = @gflops bm_func($mass_panel,$U_panel,$V_panel)



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

# mass_panel_ref(u,v) =  ∫( ( u*v )*meas_ref_cf )dΩ_ref
mass_panel_ref(u,v) =  biform(meas_ref_cf,dΩ_ref)(u,v)

# @benchmark  bm_func($mass_panel_ref,$U_panel,$V_panel)
# t_ref = @belapsed bm_func($mass_panel_ref,$U_panel,$V_panel)
ops_ref = @count_ops bm_func($mass_panel_ref,$U_panel,$V_panel)
flops_ref = @gflops bm_func($mass_panel_ref,$U_panel,$V_panel)

total_counts(ops_ambient)
total_counts(ops_panel)
total_counts(ops_ref)

0.83 GFlops,  0.45% peak  (8.59e+05 flop, 1.03e-03 s, 18793 alloc: 1.16 MiB)
2.28 GFlops,  1.16% peak  (7.12e+05 flop, 3.13e-04 s, 19871 alloc: 429.84 KiB)
2.29 GFlops,  1.20% peak  (6.81e+05 flop, 2.98e-04 s, 18949 alloc: 407.97 KiB)

0.83e8*1.03e-03
2.28e8*3.13e-04
2.29e8 *2.98e-04
