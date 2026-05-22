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
ambient_model = SingleAmbientDiscreteModel(panel_model, 1.0)

model = ambient_model
trian = Triangulation(model)

### Integrate on the physical domain
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

### standard grad grad term
Iₖ = lazy_map(Broadcasting(Operation(⋅)),∇ϕₖ,∇ϕₖᵀ)


Qₕ_cell_point = get_cell_points(Qₕ)
qₖ = get_data(Qₕ_cell_point)

J = lazy_map(∇,cmaps)
Jq = lazy_map(evaluate,J,qₖ)
intq = lazy_map(evaluate,Iₖ,qₖ)
iwq = lazy_map(IntegrationMap(),intq,Qₕ.cell_weight,Jq)

arr = iwq
cache = array_cache(arr);
ops_ambient = @count_ops lazy_collect($cache,$arr)
total_counts(ops_ambient)

arr = iwq
cache = array_cache(arr);
@btime lazy_collect($cache,$arr)

arr = iwq
cache = array_cache(arr);
@gflops lazy_collect($cache,$arr)
# 0.11 GFlops,  0.15% peak  (1.00e+00 flop, 9.38e-09 s, 0 alloc: 0 bytes)

# 3961


### make sure all array caches evaluate
function Gridap.Arrays.array_cache(dict::Dict,a::Gridap.Arrays.LazyArray)
  println("here")
  Gridap.Arrays._array_cache!(dict,a)
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

Iₖ = Broadcasting(Operation(⋅))(∇ϕₖ,∇ϕₖᵀ)

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

ops_ambient = @count_ops bm_intergrate($intq,$w,$Jq)
total_counts(ops_ambient)
932

@gflops bm_intergrate($intq,$w,$Jq)
# 0.73 GFlops,  0.99% peak  (9.32e+02 flop, 1.28e-06 s, 28 alloc: 2.38 KiB)


"""
redo ...
"""

################################################################################
∇ϕr = map(x->∇(x),ϕr )
∇ϕrₖ = ∇ϕr
manual_grad_dv_array = Broadcasting(push_∇)(∇ϕrₖ,cmaps[1])

∇ϕrᵀ                 = map(x->∇(transpose(x)),(ϕr))
∇ϕrₖᵀ                = ∇ϕrᵀ
manual_grad_du_array = Broadcasting(push_∇)(∇ϕrₖᵀ,cmaps[1])


∇ϕₖ = manual_grad_dv_array
∇ϕₖᵀ = manual_grad_du_array

Iₖ = Broadcasting(Operation(⋅))(∇ϕₖ,∇ϕₖᵀ)

qₖ = collect(get_data(Qₕ_cell_point))[1]

J = ∇(cmaps[1])
Jq = evaluate(J,qₖ)
intq = evaluate(Iₖ,qₖ)
w = Qₕ.cell_weight[1]

function bm_intergrate(intq,w,J)
  s = 0.0
  for i in 1:size(intq,2)
    s += evaluate!(cache, IntegrationMap(),intq[:,i],w,J)
  end
  return s
end

ops_ambient = @count_ops bm_intergrate($intq,$w,$Jq)
total_counts(ops_ambient)
932

@gflops bm_intergrate($intq,$w,$Jq)
# 0.73 GFlops,  0.99% peak  (9.32e+02 flop, 1.28e-06 s, 28 alloc: 2.38 KiB)
