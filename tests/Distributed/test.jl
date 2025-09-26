using DrWatson
using Gridap
using GridapDistributed
using PartitionedArrays
using MPIPreferences

using Gridap.Geometry
using Gridap.Adaptivity
MPIPreferences.use_jll_binary()

### This is Jordi's function from GridapDistributed/test/AdaptivityTests.jl
function DistributedAdaptivityGlue(serial_glue,parent,child)
  glue = map(partition(get_cell_gids(parent)),partition(get_cell_gids(child))) do parent_gids, child_gids
    old_g2l = global_to_local(parent_gids)
    old_l2g = local_to_global(parent_gids)
    new_l2g = local_to_global(child_gids)

    n2o_cell_map  = lazy_map(Reindex(old_g2l),serial_glue.n2o_faces_map[3][new_l2g])
    n2o_faces_map = [Int64[],Int64[],collect(n2o_cell_map)]
    n2o_cell_to_child_id = serial_glue.n2o_cell_to_child_id[new_l2g]
    rrules = serial_glue.refinement_rules[old_l2g]
    Gridap.Adaptivity.AdaptivityGlue(n2o_faces_map,n2o_cell_to_child_id,rrules)
  end
  return glue
end


#### The distributed panel ids are extracted from the serial. This includes both
#### owned+ghost panel_ids.
#### The owned panel_ids are extracted by determining the owned cell
function distributed_panel_ids(dmodel,spanel_ids::AbstractArray{Int})
  gids = get_cell_gids(dmodel)

  dpanel_ids = map(partition(gids)) do ids
    lid_to_gid = local_to_global(ids)
    return spanel_ids[lid_to_gid]
  end

  owned_panel_ids = map(dpanel_ids,partition(gids)) do panel_ids, cids
    owned_cells = own_to_local(cids)
    return panel_ids[owned_cells]
  end

  return dpanel_ids, owned_panel_ids
end

################################################################################
using GridapGeosciences
s_model_coarse = coarse_parametric_model()
s_model_ref = refine(s_model_coarse)
s_model_ref_ref = refine(s_model_ref)
s_model_ref_ref_ref = refine(s_model_ref_ref)


spanel_ids = [get_panel_ids(s_model_ref_ref_ref),get_panel_ids(s_model_ref_ref),get_panel_ids(s_model_ref),get_panel_ids(s_model_coarse)]

# s_model_coarse = Geometry.UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(2,3)))
# s_model_ref = refine(s_model_coarse)
# s_model_ref_ref = refine(s_model_ref)

################################################################################
##### Distributed models
################################################################################
nprocs = 6
ranks  = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end


part_to_cells = [PartitionedArrays.local_range(rank,nprocs,num_cells(s_model_coarse)) for rank in 1:nprocs]
coarse_cell_to_part = zeros(Int32,num_cells(s_model_coarse))
for (rank, cells) in enumerate(part_to_cells)
  coarse_cell_to_part[cells] .= rank
end



models = [s_model_ref_ref_ref.model, s_model_ref_ref.model, s_model_ref.model, s_model_coarse]
glues = [s_model_ref_ref_ref.glue, s_model_ref_ref.glue, s_model_ref.glue]
cell_to_part = Vector{Any}(undef,length(models))
cell_to_part[end] = coarse_cell_to_part
for level in length(models)-1:-1:1
  n2o_cells = glues[level].n2o_faces_map[3]
  cell_to_part[level] = cell_to_part[level+1][n2o_cells]
end

## plot the partition in serial
for (level,(model,cparts,panel_ids)) in enumerate(zip(models,cell_to_part,spanel_ids))
  geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
  writevtk(Triangulation(model),dir*"/serial_model_ref_$(level)";append=false,celldata=["part" => cparts],geo_map=geo_map)
end

dmodels = Vector{Any}(undef,length(models))
dpanel_ids = Vector{Any}(undef,length(models))
owned_panel_ids = Vector{Any}(undef,length(models))

dmodels[end] = DiscreteModel(ranks,models[end],cell_to_part[end])
dpanel_ids[end],owned_panel_ids[end] = distributed_panel_ids(dmodels[end],spanel_ids[end])
for level in length(models)-1:-1:1
  child = DiscreteModel(ranks,models[level],cell_to_part[level])
  parent = dmodels[level+1]
  glue = DistributedAdaptivityGlue(glues[level],parent,child)
  dmodels[level] = GridapDistributed.DistributedAdaptedDiscreteModel(child,parent,glue)
  dpanel_ids[level],owned_panel_ids[level] = distributed_panel_ids(child,spanel_ids[level])
end

for i in 1:length(models)
  println("Ref lvl: ", length(models)-i)
  map(local_views(dpanel_ids[i]),local_views(owned_panel_ids[i]),ranks) do p,op, r
    println("Proc no: $r")
    println("\tOwned+Ghost panel_ids: ", p)
    println("\tOwned panel_ids: ", op)
  end
end

#### Plot models 1 at a time. If using a loop, vtk crashes
include("vtk.jl")

level = 1
dmodel = dmodels[level]
o_pids = owned_panel_ids[level]

cell_geo_map = map(o_pids) do pid
  return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pid)
end

writevtk(Triangulation(dmodel),dir*"/ambient_model_ref_$(level)", append=false, compress=false,geo_map=cell_geo_map)



map(local_views(dmodel),ranks,o_pids) do model, r, pid
  println(typeof(model.model))
  grid = get_grid(model)

  cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pid)
  writevtk(Triangulation(model.model),dir*"/ambient_model_ref_1_rank$r", append=false, compress=false,geo_map=cell_geo_map)

end


##################################################################################
###### Laplacian test
################################################################################
level = 1
panel_model = dmodels[level]
panel_ids = owned_panel_ids[level]
p_fe = 1

function f(p)
  function _f(αβ)
    α,β = αβ
    if p == 1 || p == 5 || p == 6
      return RADIUS^3/rho3(αβ)*tan(α)*tan(β)
    else
      return -RADIUS^3/rho3(αβ)*tan(α)*tan(β)
    end
  end
end


Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,2*p_fe+1)

V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1, constraint=:zeromean)
U = TrialFESpace(V)


function distributed_panelwise_cellfield(f::Function,
    trian::GridapDistributed.DistributedTriangulation,
    panel_ids::AbstractArray)
  fields = map(trian.trians,panel_ids) do trian, pids
      panelwise_cellfield(f,trian,pids)
  end
  GridapDistributed.DistributedCellField(fields,trian)
end

f_panel_cf = distributed_panelwise_cellfield(f,Ω_panel,panel_ids)

function geo_map(panel_ids::AbstractArray{Int})
  return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
end

function distributed_geo_map(owned_panel_ids::AbstractArray{Int})
   cell_geo_map = map(owned_panel_ids) do pid
    return geo_map(pid)
  end
  return cell_geo_map
end

cell_geo_map = map(o_pids) do pid
  return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pid)
end

writevtk(Triangulation(panel_model),dir*"/ambient_model_cf_ref_$(level)",cellfields=["f"=>f_panel_cf],
append=false,geo_map=cell_geo_map)



inv_metric_cf = CellField(analytic_inv_metric,Ω_panel)
meas_cf = CellField(sqrtg,Ω_panel)
slap_panel_cf =  distributed_panelwise_cellfield(surflap(f),Ω_panel,panel_ids)

println("Zeromean: ", sum(∫(f_panel_cf*meas_cf)dΩ))
@check sum(∫(f_panel_cf*meas_cf)dΩ) < 1e-14 "Function must be zero mean to solve with zeromean FE space!"

rhs_cf = - slap_panel_cf

poisson_biform(u,v) =  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
poisson_liform(v) = ∫(  (rhs_cf*v)*meas_cf )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

uh = solve(LUSolver(),op)

e = l2(f_panel_cf-uh,dΩ)
eh = f_panel_cf-uh

writevtk(Ω_panel,dir*"/ambient_model_nref$(level)_p$p_fe",cellfields=["eu"=>eh],append=false,geo_map=cell_geo_map)
