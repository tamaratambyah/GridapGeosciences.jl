using Gridap
using GridapGeosciences
using DrWatson
using FillArrays
using Gridap.Geometry, Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Test


using GridapDistributed
using PartitionedArrays
using MPIPreferences
MPIPreferences.use_jll_binary()


dir = datadir("Distributed")
include("../convergence_tools.jl")


n_ref_lvls = 2

nprocs = 6

ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end


s_models  = get_refined_models(n_ref_lvls,true)

dmodels, = get_distributed_refined_models(ranks,nprocs,s_models)

################################################################################
#### Serial models
################################################################################
model = s_models[2]
trian = Triangulation(model)
get_panel_ids(trian)

trian = BoundaryTriangulation(model)
get_panel_ids(trian)

trian = SkeletonTriangulation(model)
get_panel_ids(trian)

################################################################################
#### BodyFittedTriangulation
################################################################################

include("../Advection/advection_funcs.jl")
u = panel_to_cartesian(u0)
vX = panel_to_cartesian(tangent_vec(vecX))
uvX = panel_to_cartesian(u0vecX)

dmodel = dmodels[2]
get_panel_ids(dmodel)


gids = get_cell_gids(dmodel)
trian = Triangulation(dmodel)
get_panel_ids(trian)


panelwise_cellfield(u,trian,get_panel_ids(trian))

map(trian.trians,partition(gids)) do t, cid
  panel_ids = get_panel_ids(t)
  owned_cells = own_to_local(cid)
  panel_ids[owned_cells]
end

Ω_panel = Triangulation(dmodel)
panel_ids = get_panel_ids(dmodel)

_rhs(p) = αβ -> u(p)(αβ) + surfdiv(contra_v(uvX))(p)(αβ)

v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)
rhs_cf = panelwise_cellfield(_rhs,Ω_panel,panel_ids)

p_fe = 1
Q = TestFESpace(dmodel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)

# hard code RT space as order 1 -- for velocity
V = TestFESpace(dmodel, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
U = TrialFESpace(V)

uh = interpolate(u_cf,P)

cell_geo_map = geo_map_func(get_owned_panel_ids(dmodel))
labels = ["uh","rhs","v_contr_cf"]
panel_cfs = [uh,rhs_cf,v_contr_cf]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/distributed_model", cellfields=cellfields,append=false,geo_map=cell_geo_map)



#### boundary cf
fids = get_face_gids(dmodel,2)
btrian = BoundaryTriangulation(dmodel)
get_panel_ids(btrian)

panelwise_cellfield(u,btrian)


map(btrian.trians,partition(fids),panel_ids) do t, fid,pid
  panel_ids = get_panel_ids(t)
  owned_faces = own_to_local(fid)
  panel_ids
end

strian = SkeletonTriangulation(dmodel)
get_panel_ids(strian)
sids = map(strian.trians) do t
  get_panel_ids(t)
end




map(trian.trians) do trian
  cmap = get_cell_map(trian)
  pts = get_cell_ref_coordinates(trian)
  lazy_map(evaluate,cmap,pts)
  @test true
end

cell_geo_map = geo_map_func(get_owned_panel_ids(trian))
writevtk(trian,dir*"/distributed_model",append=false,geo_map=cell_geo_map)

u_cf = panelwise_cellfield(u,trian)

################################################################################
#### BoundaryTriangulation
#### Need to return the mask that is all the interior cells
################################################################################

btrian = BoundaryTriangulation(dpanel_model)
map(btrian.trians) do t
  get_panel_ids(t)
end



model = get_background_model(btrian)

fids = get_face_gids(model,2)
dpanel_ids = get_panel_ids(model)

panel_ids = map(dpanel_ids,partition(fids)) do pids, fids
  owned_cells = own_to_local(fids)
  return pids[owned_cells]
end


b_panel_ids =  get_panel_ids(btrian)

b_geo_map = geo_map_func(b_panel_ids)
writevtk(btrian,dir*"/distributed_boundary_trian",append=false,geo_map=b_geo_map)



################################################################################
#### SkeletonTriangulation
#### Need to dispatch to serial
################################################################################
dpanel_model = dmodels[2]
skel = SkeletonTriangulation(dpanel_model)

get_panel_ids(skel)

model = get_background_model(skel)
dpanel_ids = get_panel_ids(model)
gids = get_face_gids(model,2)
map(partition(gids),dpanel_ids) do g,pids
#  return pids[own_to_local(g)]
println(own_to_local(g))
end

proc = 2
ttrian = skel.trians.items[proc].plus

trian = ttrian.parent
_panel_ids = get_panel_ids(model)
panel_ids = _panel_ids.items[proc]
Dc = 2
glue = get_glue(trian,Val(Dc))
face_2_cell = glue.tface_to_mface
face_panel_ids = panel_ids[face_2_cell]

Gridap.Geometry.restrict(face_panel_ids,ttrian.cell_to_parent_cell)





panel_model = get_background_model(btrian)
  panel_ids = get_panel_ids(panel_model)
  println(length(panel_ids))
  Dc = num_cell_dims(panel_model)

  glue = get_glue(btrian,Val(Dc))
  face_2_cell = glue.tface_to_mface
  face_panel_ids = panel_ids[face_2_cell]
  return face_panel_ids


skel_panel_ids = get_panel_ids(skel)
_skel_panel_ids = get_skel_panel_ids(skel_panel_ids)

skel_geo_map = geo_map_func(_skel_panel_ids)
writevtk(skel,dir*"/distributed_skel_trian",append=false,geo_map=skel_geo_map)

n_Λ = get_normal_vector(skel)
writevtk(skel,dir*"/distributed_skel_trian",cellfields=["np"=>n_Λ.plus,"nm"=>n_Λ.minus,"diffn"=>n_Λ.plus+n_Λ.minus],append=false,geo_map=skel_geo_map)

######
dpanel_model = dmodels[2]
trian = SkeletonTriangulation(dpanel_model)
pts = get_cell_points(trian)

############# area form
area_form = pullback_area_form(trian)

cf1 = area_form.plus(pts)
cf2 = area_form.minus(pts)
map(cf1,cf2) do c1,c2
  @test length(c1) == length(c2)
  @test sum(c1 .≈ c2) == length(c1)
end

#### facet normal
cell_geo_map = geo_map_func(get_panel_ids(dpanel_model))
function Geometry.get_facet_normal(trian::Gridap.Geometry.TriangulationView,cell_geo_map::AbstractArray)
  println("trian views")
   Geometry.get_facet_normal(trian.parent,cell_geo_map)
end

fields = map(trian.trians,cell_geo_map) do t, m
  return get_facet_normal(t,m)
end
GridapDistributed.DistributedCellField(fields,trian)
typeof(fields) <:AbstractArray{<:SkeletonPair}






##### push forward normal
function GridapGeosciences.pushforward_normal(trian::Gridap.Geometry.TriangulationView)
  _pushforward_normal(trian)
end


fields = map(trian.trians) do t
  n_mapped = pushforward_normal(t)
  return n_mapped
end
n_mapped = GridapDistributed.DistributedCellField(fields,trian)

pts = get_cell_points(trian)
(n_mapped.plus - n_mapped.minus)(pts)
@test sum(n_mapped.plus(pts) .≈ n_3D.plus(pts)) == num_facets(panel_model)



######################
