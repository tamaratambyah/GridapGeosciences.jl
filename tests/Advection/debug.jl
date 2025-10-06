using Gridap
using GridapGeosciences

include("../convergence_tools.jl")
include("advection_funcs.jl")

using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Gridap.Helpers

n_ref_lvls = 2

vX = panel_to_cartesian(tangent_vec(vecX))

models  = get_refined_models(n_ref_lvls,true)

using Test
using GridapDistributed
using PartitionedArrays
using MPIPreferences
MPIPreferences.use_jll_binary()

nprocs = 6

ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end

dmodels,  = get_distributed_refined_models(ranks,nprocs,models)

panel_model = dmodels[1]
Ω_panel = Triangulation(panel_model)

trian = Ω_panel



trian_panel_ids = get_panel_ids(trian)
model_panel_ids = get_panel_ids(panel_model)
map(trian_panel_ids,model_panel_ids) do t, m
  @test t == m
end





for panel_model in models
  Ω_panel = Triangulation(panel_model)
  @test get_panel_ids(Ω_panel) == get_panel_ids(panel_model)
end

f = contra_v(vX)
panel_ids = get_panel_ids(trian)
fields = map(trian.trians,panel_ids) do trian,pid
  cell_field = map(p->GenericField(f(p)),pid)
  Gridap.CellData.GenericCellField(cell_field,trian,PhysicalDomain())
end
cf = GridapDistributed.DistributedCellField(fields,trian)
vel = interpolate(cf,U)

# func(αβ) = VectorValue(αβ[1],αβ[2])
# cf = CellField(func,trian)
# interpolate(cf,U)

v_contr_cf = panelwise_cellfield(contra_v(vX),Ω_panel)


Λ = SkeletonTriangulation(panel_model)

Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)

# hard code RT space as order 1 -- for velocity
V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
U = TrialFESpace(V)


trian = Ω_panel
f = contra_v(vX)
fields = map(trian.trians,panel_ids) do trian, pids
  println(typeof(trian.trian)<:Gridap.Geometry.BodyFittedTriangulation)
  cell_field = map(p->Gridap.Fields.GenericField(f(p)),pids)
  Gridap.CellData.GenericCellField(cell_field,trian.trian,PhysicalDomain())
end
v_contr_cf = GridapDistributed.DistributedCellField(fields,trian)
vel = interpolate(v_contr_cf,U)


_a(u,v) = ∫( u⋅v )dΩ
_l(v) = ∫( v_contr_cf ⋅ v )dΩ
op = AffineFEOperator(_a,_l,U,V)
vel = solve(LUSolver(),op)
