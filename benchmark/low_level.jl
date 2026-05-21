"""
In this test, benchmark the assembly of the grad-grad term in the parametric model
vs. the ambient model on a single panel
"""

using GridapGeosciences
using Gridap
using BenchmarkTools
using DrWatson
using Gridap.Geometry, Gridap.CellData, Gridap.Arrays, Gridap.FESpaces
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



nC = 2^nref
_panel_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(nC,nC))
panel_model = SingleParametricDiscreteModel(_panel_model)
ambient_model = SingleAmbientDiscreteModel(panel_model,radius)

degree = 5*(p_fe+1)

################################################################################
########## Ambient model
################################################################################
Ω_ambient = Triangulation(ambient_model)
dΩ_ambient = Measure(Ω_ambient,degree)
U_ambient, V_ambient = set_up(ambient_model,reffee,conf)


function collect_measure(meas)
  quad = meas.quad
  CellQuadrature(
    collect(quad.cell_quad),
    collect(quad.cell_point),
    collect(quad.cell_weight),
    quad.trian,
    quad.data_domain_style,
    quad.integration_domain_style
  ) |> Measure
end

function collect_cellfield(u)
  GenericCellField(collect(get_data(u)),get_triangulation(u),DomainStyle(u))
end

function Gridap.Arrays.array_cache(dict::Dict,a::Gridap.Arrays.LazyArray)
  Gridap.Arrays._array_cache!(dict,a)
end

function lazy_collect(cache,arr)
  s = 0.0
  for i in eachindex(arr)
    ai = getindex!(cache,arr,i)
    s += ai[1,1]
  end
  return s
end

u = get_trial_fe_basis(U_ambient)
v = get_fe_basis(V_ambient)
dΩ_ambient_2 = collect_measure(dΩ_ambient)
contr = ∫(∇(u)⋅ ∇(v))dΩ_ambient_2

arr = get_array(contr);
open(joinpath(pwd(),"ambient_tree.txt"),"w") do f
  print_op_tree(f,arr;maxdepth=20)
end
cache = array_cache(arr);

@gflops lazy_collect($cache,$arr)
ops = @count_ops lazy_collect($cache,$arr)




################################################################################
########## Parametric model -- integrate in the physical domain (includes jacobians)
################################################################################
Ω_panel = Triangulation(panel_model)
dΩ_panel = Measure(Ω_panel,degree)
dΩ_panel_2 = collect_measure(dΩ_panel)
U_panel, V_panel = set_up(panel_model,reffee,conf)

u = get_trial_fe_basis(U_panel)
v = get_fe_basis(V_panel)
dΩ_panel_2 = collect_measure(dΩ_panel)

# inv_metric = CellField(inv_g,Ω_panel)
# meas_cf = CellField(_sqrtg,Ω_panel)

# contr = ∫(  ( ∇(v)⋅ (inv_metric⋅ ∇(u) ) )*meas_cf)dΩ_panel_2

allg_cf = CellField(allg,Ω_panel)
contr = ∫(  ( ∇(v)⋅ (allg_cf⋅ ∇(u) ) ))dΩ_panel_2

arr = get_array(contr)
open(joinpath(pwd(),"ref_tree.txt"),"w") do f
  print_op_tree(f,arr;maxdepth=20)
end
cache = array_cache(arr)

ops = @gflops lazy_collect($cache,$arr)
ops = @count_ops lazy_collect($cache,$arr)


################################################################################
########## Parametric model -- integrate in reference domain (no jacobians)
########## Now the metric goes from the refee -> ambient space
################################################################################
quad = CellQuadrature(Ω_panel,degree;
        data_domain_style=ReferenceDomain(),
        integration_domain_style=ReferenceDomain())
dΩ_ref = Gridap.CellData.GenericMeasure(quad)
dΩ_ref_2 = collect_measure(dΩ_ref)

cmap = get_cell_map(panel_model.grid)
g_ref = lazy_map(Operation(allg),cmap)
g_ref_cf = Gridap.CellData.GenericCellField(g_ref,Ω_panel,ReferenceDomain())


contr = ∫(  ( ∇(v)⋅ (g_ref_cf⋅ ∇(u) ) ))dΩ_ref_2

arr = get_array(contr)
cache = array_cache(arr)
ops = @gflops lazy_collect($cache,$arr)
ops = @count_ops lazy_collect($cache,$arr)

total_counts(ops)

0.90 GFlops,  0.63% peak  (1.61e+06 flop, 1.78e-03 s, 34976 alloc: 1.62 MiB)
0.77 GFlops,  0.42% peak  (3.39e+05 flop, 4.41e-04 s, 11008 alloc: 227.75 KiB)
0.74 GFlops,  0.41% peak  (3.39e+05 flop, 4.58e-04 s, 10976 alloc: 226.00 KiB)


1399856
310800
310800
;
