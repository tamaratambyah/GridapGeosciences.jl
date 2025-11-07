using DrWatson
using Gridap
using GridapDistributed
using GridapP4est
using GridapGeosciences
using GridapGeosciences.Distributed

using MPI
using PartitionedArrays


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

num_horizontal_uniform_refinements = 2
num_vertical_uniform_refinements = 1
omodel = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
	                                       num_horizontal_uniform_refinements=num_horizontal_uniform_refinements,
                                           num_vertical_uniform_refinements=num_vertical_uniform_refinements);


dir = datadir("Distributed")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

include("../convergence_tools.jl")
include("missing_overloads.jl")
include("williamson_funcs_3D.jl")

vX = panel_to_cartesian(tangent_vec(u₀(0.0)))
p_fe = 1
panel_model = omodel.parametric_dmodel
panel_ids = get_panel_ids(panel_model)
tags = ["bottom_boundary",  "top_boundary"]

###### no ghost
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,4)
V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv, dirichlet_tags=tags)
U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)

covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
u_contra_h = interpolate(u_contra_cf,U)

u_proj_h = covarient_basis_cf ⋅ u_contra_cf
e = l2(u_proj_cf-u_proj_h,dΩ)
i_am_main(ranks) && println(e)

cell_geo_map = geo_map_func(Ω_panel)
panel_cfs = [u_proj_cf, u_proj_h,u_proj_h-u_proj_cf]
labels = ["u", "uh","eu"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/no_ghost",cellfields=cellfields,append=false,geo_map=cell_geo_map)

### with ghosts
Ω_panel = Triangulation(with_ghost,panel_model)
dΩ = Measure(Ω_panel,4)
V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv, dirichlet_tags=tags)
# b_cf = CellField(VectorValue(0.0,0.0,0.0),Ω_panel)
U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)

covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
u_contra_h = interpolate(u_contra_cf,U)

u_proj_h = covarient_basis_cf ⋅ u_contra_cf

e = l2(u_proj_cf-u_proj_h,dΩ)
i_am_main(ranks) && println(e)
cell_geo_map = geo_map_func(get_panel_ids(Ω_panel))
panel_cfs = [u_proj_cf, u_proj_h,u_proj_h-u_proj_cf]
labels = ["u", "uh","eu"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/with_ghost",cellfields=cellfields,append=false,geo_map=cell_geo_map)
