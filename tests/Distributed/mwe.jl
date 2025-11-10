using DrWatson
using Gridap
using GridapDistributed
using GridapP4est
using GridapGeosciences

using MPI
using PartitionedArrays


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

include("../Geophysical/Williamson2Test.jl")

dir = datadir("DistributedWaveEquation")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

h = panel_to_cartesian(h₀(0.0))
vX = panel_to_cartesian(tangent_vec(u₀(0.0)))
p_fe = 1
ls = LUSolver()

omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=3)
panel_model = omodel.parametric_dmodel

# das = FullyAssembledRows()
das = SubAssembledRows()

panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(das,panel_model)
dΩ = Measure(Ω_panel,2*(p_fe+1))

### See https://github.com/tamaratambyah/GridapGeosciences.jl/issues/5
### use triangulation for FE space
Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)

V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
U = TrialFESpace(V)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)
u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

# interpolate into FE space
xh_interpolate = interpolate([u_contra_cf,h_cf],X)
u_contra_h = xh_interpolate[1] ## contravariant components
u_interpolate = covarient_basis_cf⋅u_contra_h ## physical velocity
h_interpolate = xh_interpolate[2]

sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(vX)),Ω_panel,panel_ids)
sgrad_cf = panelwise_cellfield(sgrad(h),Ω_panel,panel_ids)
pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω_panel,panel_ids)

# manufacture rhs functions
rhs_scalar = h_cf + sdiv_cf
rhs_vector = u_proj_cf + sgrad_cf
rhs_con_vector = pinvJ_cf ⋅ rhs_vector # exact contravariant component


metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)

biform1((u,p),(v,q)) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
biform2((u,p),(v,q)) = ∫( (p*q)*meas_cf )dΩ + ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ

biformX((u,p),(v,q)) = biform1((u,p),(v,q)) + biform2((u,p),(v,q))
liformX((v,q)) = ∫( rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ + ∫( (rhs_scalar*q)*meas_cf )dΩ

assem = SparseMatrixAssembler(X,Y,das)
op = AffineFEOperator(biformX,liformX,X,Y,assem)
uh,ph = solve(ls,op)

uh_proj = covarient_basis_cf ⋅ uh

Ω_error = Triangulation(panel_model)
dΩ_error = Measure(Ω_error,6*p_fe+1)
e_u = l2( (u_proj_cf - uh_proj),meas_cf,dΩ_error) # error in physical velocity u
e_p = l2((h_cf - ph),meas_cf,dΩ_error) # error in depth



_Ω_panel = Triangulation(panel_model)
cell_geo_map = geo_map_func(_Ω_panel)

if das == FullyAssembledRows()
  cell_geo_map = geo_map_func(get_panel_ids(_Ω_panel))
end

panel_cfs = [h_cf, u_proj_cf, ph, uh_proj, h_cf-ph, u_proj_cf-uh_proj, u_interpolate, h_interpolate ]
labels = ["p","u_proj", "ph", "uh_proj", "ep","eu", "u_interpolate", "h_interpolate"]

cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/sol_proc$(length(ranks))",cellfields=cellfields,append=false,geo_map=cell_geo_map)
