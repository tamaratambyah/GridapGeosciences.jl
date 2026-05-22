using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Arrays, Gridap.FESpaces
using Gridap.ReferenceFEs, Gridap.Fields
using FillArrays
using GridapGeosciences
using BenchmarkTools
using GFlops

include("../single_panel_ambient_model.jl")
include("helpers.jl")



degree = 4

panel_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(1,1))


model = panel_model

trian = Triangulation(model)

## integration on the physical domain
Qₕ = CellQuadrature(trian,degree;
        data_domain_style=ReferenceDomain(),
        integration_domain_style=PhysicalDomain())


cmaps = get_cell_map(model)
reffes = get_reffes(model)
reffe = reffes[1]

ϕr = get_shapefuns(reffe)
∇ϕr = Broadcasting(∇)(ϕr)
∇ϕrₖ = Fill(∇ϕr,1)
manual_grad_dv_array = lazy_map(Broadcasting(push_∇),∇ϕrₖ,cmaps)

∇ϕrᵀ                 = Broadcasting(∇)(transpose(ϕr))
∇ϕrₖᵀ                = Fill(∇ϕrᵀ,1)
manual_grad_du_array = lazy_map(Broadcasting(push_∇),∇ϕrₖᵀ,cmaps)



∇ϕₖ = manual_grad_dv_array
∇ϕₖᵀ = manual_grad_du_array

### get the metric
allg_cf = CellField(allg,trian)
allg_array = get_data(allg_cf)
_Iₖ = lazy_map(Broadcasting(Operation(⋅)),allg_array,∇ϕₖᵀ)
Iₖ = lazy_map(Broadcasting(Operation(⋅)),∇ϕₖ,_Iₖ)


Qₕ_cell_point = get_cell_points(Qₕ)
qₖ = get_data(Qₕ_cell_point)

J = lazy_map(∇,cmaps)
Jq = lazy_map(evaluate,J,qₖ)
intq = lazy_map(evaluate,Iₖ,qₖ)
iwq = lazy_map(IntegrationMap(),intq,Qₕ.cell_weight,Jq)

arr = iwq
cache = array_cache(arr);

# @gflops lazy_collect($cache,$arr)

arr = iwq
cache = array_cache(arr);
ops_panel = @count_ops lazy_collect($cache,$arr)
total_counts(ops_panel)

arr = iwq
cache = array_cache(arr);
@btime lazy_collect($cache,$arr)


2506

## make sure all array caches evaluate
function Gridap.Arrays.array_cache(dict::Dict,a::Gridap.Arrays.LazyArray)
  println("here")
  Gridap.Arrays._array_cache!(dict,a)

  # cache = _get_cache(dict,a)
  # if cache === nothing
  #   _cache = _array_cache!(dict,a)
  #   dict[objectid(a)] = (a,_cache)
  # else
  #   _cache = cache
  # end
  # _cache
end

################################################################################



∇ϕr = map(x->∇(x),ϕr )#Broadcasting(∇)(ϕr)
∇ϕrₖ = [∇ϕr]
manual_grad_dv_array = lazy_map(Broadcasting(push_∇),∇ϕrₖ,cmaps)

∇ϕrᵀ                 = map(x->∇(transpose(x)),(ϕr))
∇ϕrₖᵀ                = [∇ϕrᵀ]
manual_grad_du_array = lazy_map(Broadcasting(push_∇),∇ϕrₖᵀ,cmaps)


∇ϕₖ = manual_grad_dv_array[1]
∇ϕₖᵀ = manual_grad_du_array[1]

allg_cf = CellField(allg,trian)
allg_array = get_data(allg_cf)[1]
_Iₖ = Broadcasting(Operation(⋅))(allg_array,∇ϕₖᵀ)
Iₖ = Broadcasting(Operation(⋅))(∇ϕₖ,_Iₖ)

qₖ = collect(get_data(Qₕ_cell_point))[1]

J = ∇(cmaps[1])
Jq = evaluate(J,qₖ)
intq = evaluate(Iₖ,qₖ)
w = Qₕ.cell_weight[1]

function bm_intergrate(intq,w,J)
  s = 0.0
  for i in 1:size(intq,2)
    s += evaluate(IntegrationMap(),intq[:,i],w,J)
  end
  return s
end

@gflops bm_intergrate($intq,$w,$Jq)
# 0.29 GFlops,  0.42% peak  (3.60e+02 flop, 1.23e-06 s, 28 alloc: 2.38 KiB)

ops_panel = @count_ops bm_intergrate($intq,$w,$Jq)
total_counts(ops_panel)
360
