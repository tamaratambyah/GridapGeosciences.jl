using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson
using Test

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

dir = datadir("Omodel_2D")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

function adapt_model(ranks,model::ParametricOctreeDistributedDiscreteModel)
  cell_partition=get_cell_gids(model.octree_dmodel)
  ref_flags=map(ranks,partition(cell_partition)) do rank,indices
      flags=zeros(Cint,length(indices))
      flags.=refine_flag
  end
  Gridap.Adaptivity.adapt(model,ref_flags)
end

include("../Laplace/analytic_funcs.jl")
include("../helpers.jl")
p_fe = 2
f = panel_to_cartesian(fX)
tol = 1e-10
ls = LUSolver()

#### coarse model
omodelH = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=1)
coarse_model = omodelH.parametric_dmodel
panel_idsH = get_panel_ids(coarse_model)
ΩH = Triangulation(coarse_model)
dΩH = Measure(ΩH,4*p_fe+1)
dΩH_error = Measure(ΩH,8*p_fe+1)

VH = TestFESpace(coarse_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
UH = TrialFESpace(VH)

f_cfH = panelwise_cellfield(f,ΩH,panel_idsH)
meas_cfH = panelwise_cellfield(sqrtg,ΩH,panel_idsH)

uH = interpolate(f_cfH,UH)
e = uH-f_cfH
el2 = l2(e,dΩH_error)
el2 = l2(e,meas_cfH,dΩH_error)


### fine model
omodelh = adapt_model(ranks,omodelH)
fine_model = omodelh.parametric_dmodel
panel_idsh = get_panel_ids(fine_model)
Ωh = Triangulation(fine_model)
dΩh = Measure(Ωh,4*p_fe+1)
dΩh_error = Measure(Ωh,8*p_fe+1)

Vh = TestFESpace(fine_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
Uh = TrialFESpace(Vh)

f_cfh = panelwise_cellfield(f,Ωh,panel_idsh)
meas_cfh = panelwise_cellfield(sqrtg,Ωh,panel_idsh)

uh = interpolate(f_cfh,Uh)
e = uh-f_cfh
el2 = l2(e,dΩh_error)
el2 = l2(e,meas_cfh,dΩh_error)


################################################################################
### Prolongation: coarse FEFunction -> fine FEFunction
################################################################################
# prolongation via interpolation:
uHh = interpolate(uH,Uh)
e = uh - uHh
# e = f_cfh - uHh
el2 = l2(e,dΩh_error)
el2 = l2(e,meas_cfh,dΩh_error)
@test el2 < tol

cell_geo_map = geo_map_func(Ωh)
writevtk(Ωh,dir*"/prolongation",cellfields=["u"=>uHh, "eu"=>uh - uHh],append=false,geo_map=cell_geo_map)


# prolongation via L2-projection
ahp(u,v) = ∫(v⋅u)*dΩh
lhp(v) = ∫(v⋅uH)*dΩh
oph = AffineFEOperator(ahp,lhp,Uh,Vh)
uHh = solve(ls,oph)
e = uh - uHh
# e = f_cfh - uHh
el2 = l2(e,dΩh_error)
el2 = l2(e,meas_cfh,dΩh_error)
@test el2 < tol

################################################################################
### Restriction: fine FEFunction -> coarse FEFunction
################################################################################
# restriction via interpolation
uhH = interpolate(uh,UH)
e = uH - uhH
# e = f_cfH - uhH
el2 = l2(e,dΩh_error)
el2 = l2(e,meas_cfh,dΩh_error)
@test el2 < tol

# restriction via L2-projection
dΩhH = Measure(ΩH,Ωh,4*p_fe+1)
aHp(u,v) = ∫(v⋅u)*dΩH_error
lHp(v) = ∫(v⋅uh)*dΩhH
oph = AffineFEOperator(aHp,lHp,UH,VH)
uhH = solve(ls,oph)
e = uH - uhH
# e = f_cfH - uhH
el2 = l2(e,dΩH_error)
el2 = l2(e,meas_cfH,dΩH_error)
@test el2 < tol

cell_geo_map = geo_map_func(ΩH)
writevtk(ΩH,dir*"/restriction",cellfields=["uhH"=>uhH, "eu"=>uH - uhH],append=false,geo_map=cell_geo_map)
