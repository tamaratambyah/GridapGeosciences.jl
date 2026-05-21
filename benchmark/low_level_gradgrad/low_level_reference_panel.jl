using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Arrays, Gridap.FESpaces
using Gridap.ReferenceFEs, Gridap.Fields
using FillArrays
using GridapGeosciences

using GFlops

include("../single_panel_ambient_model.jl")
include("helpers.jl")


degree = 4

model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(1,1))
trian = Triangulation(model)

### Integration in the reference space
Qₕ = CellQuadrature(trian,degree;
        data_domain_style=ReferenceDomain(),
        integration_domain_style=ReferenceDomain())



cmaps = get_cell_map(model)
reffes = get_reffes(model)
reffe = reffes[1]

ϕr = get_shapefuns(reffe)
∇ϕr = Broadcasting(∇)(ϕr)
∇ϕrₖ = Fill(∇ϕr,1)
# manual_grad_dv_array = lazy_map(Broadcasting(push_∇),∇ϕrₖ,cmaps)
manual_grad_dv_array = ∇ϕrₖ

∇ϕrᵀ                 = Broadcasting(∇)(transpose(ϕr))
∇ϕrₖᵀ                = Fill(∇ϕrᵀ,1)
# manual_grad_du_array = lazy_map(Broadcasting(push_∇),∇ϕrₖᵀ,cmaps)
manual_grad_du_array = ∇ϕrₖᵀ


∇ϕₖ = manual_grad_dv_array
∇ϕₖᵀ = manual_grad_du_array

# g_ref = lazy_map(Broadcasting(Operation(allg)),cmaps)

g_ref = lazy_map(Broadcasting(∘),Fill(GenericField(allg),1),cmaps)
allg_array = g_ref

_Iₖ = lazy_map(Broadcasting(Operation(⋅)),allg_array,∇ϕₖᵀ)
Iₖ = lazy_map(Broadcasting(Operation(⋅)),∇ϕₖ,_Iₖ)


Qₕ_cell_point = get_cell_points(Qₕ)
qₖ = get_data(Qₕ_cell_point)

cell_map1 = Fill(GenericField(identity),1)
cell_Jt = lazy_map(∇,cell_map1)
Jq = lazy_map(evaluate,cell_Jt,qₖ)
Gridap.TensorValues.meas(Jq[1][1])

### evaluate with jacobian since it is identity
intq = lazy_map(evaluate,Iₖ,qₖ)
iwq = lazy_map(IntegrationMap(),intq,Qₕ.cell_weight)

arr = iwq
cache = array_cache(arr);
ops_panel_ref = @count_ops lazy_collect($cache,$arr)

total_counts(ops_panel_ref)
2395

################################################################################

∇ϕr = map(x->∇(x),ϕr )
∇ϕrₖ = [∇ϕr]
manual_grad_dv_array = lazy_map(Broadcasting(push_∇),∇ϕrₖ,cmaps)

∇ϕrᵀ                 = map(x->∇(transpose(x)),(ϕr))
∇ϕrₖᵀ                = [∇ϕrᵀ]
manual_grad_du_array = lazy_map(Broadcasting(push_∇),∇ϕrₖᵀ,cmaps)


∇ϕₖ = manual_grad_dv_array[1]
∇ϕₖᵀ = manual_grad_du_array[1]

allg_cf = CellField(allg,trian)
g_ref = Broadcasting(Operation(allg))(cmaps[1])
allg_array = g_ref

_Iₖ = Broadcasting(Operation(⋅))(allg_array,∇ϕₖᵀ)
Iₖ = Broadcasting(Operation(⋅))(∇ϕₖ,_Iₖ)


qₖ = collect(get_data(Qₕ_cell_point))[1]

J = ∇(cmaps[1])
Jq = evaluate(J,qₖ)
intq = evaluate(Iₖ,qₖ)
w = Qₕ.cell_weight[1]

function bm_intergrate(intq,w)
  s = 0.0
  for i in 1:size(intq,2)
    s += evaluate(IntegrationMap(),intq[:,i],w)
  end
  return s
end

ops_panel_ref = @count_ops bm_intergrate($intq,$w)
total_counts(ops_panel_ref)
