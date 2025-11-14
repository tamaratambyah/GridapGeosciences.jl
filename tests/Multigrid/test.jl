using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson
using Test
import GridapDistributed: i_am_in

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
p_fe = 1
f = panel_to_cartesian(fX)
tol = 1e-6
ls = LUSolver()

#### coarse model
omodelH = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=0)
coarse_model = omodelH.parametric_dmodel
panel_idsH = get_panel_ids(coarse_model)
ΩH = Triangulation(coarse_model)
dΩH = Measure(ΩH,2*p_fe+1)

VH = TestFESpace(coarse_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
UH = TrialFESpace(VH)

f_cfH = panelwise_cellfield(f,ΩH,panel_idsH)
meas_cfH = panelwise_cellfield(sqrtg,ΩH,panel_idsH)

uH = interpolate(f_cfH,UH)
e = uH-f_cfH
el2 = l2(e,dΩH)
el2 = l2(e,meas_cfH,dΩH)


### fine model
omodelh = adapt_model(ranks,omodelH)
fine_model = omodelh.parametric_dmodel
panel_idsh = get_panel_ids(fine_model)
Ωh = Triangulation(fine_model)
dΩh = Measure(Ωh,2*p_fe+1)

Vh = TestFESpace(fine_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
Uh = TrialFESpace(Vh)

f_cfh = panelwise_cellfield(f,Ωh,panel_idsh)
meas_cfh = panelwise_cellfield(sqrtg,Ωh,panel_idsh)

uh = interpolate(f_cfh,Uh)
e = uh-f_cfh
el2 = l2(e,dΩh)
el2 = l2(e,meas_cfh,dΩh)


################################################################################
### Prolongation: coarse FEFunction -> fine FEFunction
################################################################################
# prolongation via interpolation:
uHh = interpolate(uH,Uh)
e = uh - uHh
el2 = l2(e,dΩh)
el2 = l2(e,meas_cfh,dΩh)
@test el2 < tol

# prolongation via L2-projection
ahp(u,v) = ∫(v⋅u)*dΩh
lhp(v) = ∫(v⋅uH)*dΩh
oph = AffineFEOperator(ahp,lhp,Uh,Vh)
uHh = solve(ls,oph)
e = uh - uHh
el2 = l2(e,dΩh)
el2 = l2(e,meas_cfh,dΩh)
@test el2 < tol

################################################################################
### Restriction: fine FEFunction -> coarse FEFunction
################################################################################
# restriction via interpolation
uhH = interpolate(uh,UH)
e = uH - uhH
el2 = l2(e,dΩh)
el2 = l2(e,meas_cfh,dΩh)
@test el2 < tol

# restriction via L2-projection
dΩhH = Measure(ΩH,Ωh,2*p_fe+1)
aHp(u,v) = ∫(v⋅u)*dΩH
lHp(v) = ∫(v⋅uh)*dΩhH
oph = AffineFEOperator(aHp,lHp,UH,VH)
uhH = solve(ls,oph)
e = uH - uhH
el2 = l2(e,dΩH)
el2 = l2(e,meas_cfH,dΩH)
@test el2 < tol

################################################################################
### Redistribution
################################################################################
function Gridap.Adaptivity.get_model(model::GridapDistributed.DistributedAdaptedDiscreteModel)
  println("my func ")
  GridapDistributed.GenericDistributedDiscreteModel(
    map(get_model,local_views(model)),
    get_cell_gids(model);
    metadata = nothing
  )
end
fmodel_red, red_glue = GridapDistributed.redistribute(fine_model)
# fmodel_red, red_glue = GridapDistributed.redistribute(omodelh)

# panel_idred = get_panel_ids(fmodel_red)
# Ωhred  = Triangulation(fmodel_red)
# dΩhred = Measure(Ωhred,degree)

# Vhred = TestFESpace(fmodel_red, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
# Uhred = TrialFESpace(Vhred)

# f_cfred = panelwise_cellfield(f,Ωred,panel_idsred)
# meas_cfred = panelwise_cellfield(sqrtg,Ωred,panel_idsred)

# uhred = interpolate(f_cfred,Uhred)
# e = uhred-f_cfhred
# el2 = l2(e,dΩhred)
# el2 = l2(e,meas_cfred,dΩhred)


# uhred2 = GridapDistributed.redistribute_fe_function(uh,Vhred,fmodel_red,red_glue)
# e = f_cfh - uhred2
# el2 = l2(e,dΩhred)
# el2 = l2(e,meas_cfred,dΩhred)
