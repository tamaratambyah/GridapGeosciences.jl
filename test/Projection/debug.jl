using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson

## should return ghost+owned panel ids
MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

dir = datadir("discontinuous")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)


omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=2)
panel_model = omodel.parametric_dmodel

Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,8)
panel_ids = get_panel_ids(panel_model)



function fS(p)
  function f(αβ)
    xyz = forward_map_2D(p)(αβ)
    xyz[1]*xyz[2]*xyz[3]
  end
end

vX(xyz) = VectorValue(-xyz[2], xyz[1], 0)
fV = panel_to_cartesian(tangent_vec(vX))

fS_cf = panelwise_cellfield(fS,Ω_panel,panel_ids)
fV_contra_cf = panelwise_cellfield(contra_v(fV),Ω_panel,panel_ids)

covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
fV_proj = covarient_basis_cf ⋅ fV_contra_cf

cell_geo_map = latlon_geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/model",
  cellfields=["fS"=>fS_cf,"fV"=> fV_proj, "fV_con"=>fV_contra_cf],
  append=false, geo_map=cell_geo_map)



p_fe = 1
W = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
R = TrialFESpace(W)

V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
U = TrialFESpace(V)

P = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
H = TrialFESpace(P)

L = TestFESpace(Ω_panel, ReferenceFE(lagrangian,VectorValue{2,Float64},p_fe); conformity=:H1)
K = TrialFESpace(L)

fS_L2 = interpolate(fS_cf,R)
fS_H1 = interpolate(fS_cf,H)

fV_contra_RT = interpolate(fV_contra_cf,U)
fV_contra_H1 = interpolate(fV_contra_cf,K)

metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
eu =  fV_contra_RT - fV_contra_cf
sum(∫( eu⋅(metric_cf⋅eu)*meas_cf )dΩ)

eu =  fV_contra_H1 - fV_contra_cf
sum(∫( eu⋅(metric_cf⋅eu)*meas_cf )dΩ)

writevtk(Ω_panel,dir*"/model",
  cellfields=["fS"=>fS_cf,"fV"=> fV_proj, "fV_con"=>fV_contra_cf,
              "efS_L2"=>fS_L2-fS_cf,
              "efS_H1"=>fS_H1-fS_cf,
              "efV_RT"=>covarient_basis_cf⋅fV_contra_RT -fV_proj,
              "efV_H1"=> covarient_basis_cf⋅fV_contra_H1 -fV_proj ],
  append=false, geo_map=latlon_geo_map_func(Ω_panel))
