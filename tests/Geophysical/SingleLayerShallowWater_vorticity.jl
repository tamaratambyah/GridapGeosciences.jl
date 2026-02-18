using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Test
using GridapPETSc

dir = datadir("SW_3D")
!isdir(dir) && mkdir(dir)

include("../convergence_tools.jl")
include("Williamson2Test_3D.jl")

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

n_ref_lvls = 3
ζ = 0.0

vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
f_vec = panel_to_cartesian(f_vec₀(ζ))
η_vec = panel_to_cartesian(η_vec₀(ζ))


o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
        num_horizontal_uniform_refinements=n_ref_lvls,
        num_vertical_uniform_refinements=0)
panel_model = o3model.parametric_dmodel

ls = LUSolver()
p_fe = 1

das =  FullyAssembledRows()
# das =  SubAssembledRows()

## finite element solver
panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(das,panel_model)
dΩ = Measure(Ω_panel,4*(p_fe+1))
Ω_error = Triangulation(panel_model)
dΩ_error = Measure(Ω_error,4*(p_fe+1))

Λ_panel = SkeletonTriangulation(das,panel_model)
dΛ = Measure(Λ_panel,4*(p_fe+1))


inv_jacobian(p) = x -> inv(forward_jacobian_3D(p)(x))
contra_v_3D(vecX::Function,p::Int) = αβ -> inv_jacobian(p)(αβ) ⋅ vecX(p)(αβ)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)

u_contra_cf = panelwise_cellfield(contra_v_3D(vX),Ω_panel,panel_ids)
f_cf = panelwise_cellfield(contra_v_3D(f_vec),Ω_panel,panel_ids)
η_cf = panelwise_cellfield(contra_v_3D(η_vec),Ω_panel,panel_ids)


tags = ["top_boundary", "bottom_boundary"]


V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
U = TrialFESpace(V,u_contra_cf)

R = TestFESpace(panel_model, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
H = TrialFESpace(R,η_cf)

## metric information
detg_cf = panelwise_cellfield(detg,Ω_panel,panel_ids)
metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)
covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
jac_cf = panelwise_cellfield(forward_jacobian,Ω_panel,panel_ids)





## push forward to ambient normal vector
n = CellField(VectorValue(1.0,0.0,0.0),Ω_panel)
ff = Operation(sqrt)(  n  ⋅ (inv_metric_cf⋅ n )  )
n_ambient  = (jac_cf ⋅(inv_metric_cf  ⋅n ) )/ff


writevtk(Ω_panel,dir*"/boundary.vtu", cellfields=["u"=>n×(metric_cf⋅u_contra_cf)],
      append=false,geo_map=latlon_geo_map_func(Ω_panel))


# initial conditions
u_contra_h = interpolate(u_contra_cf,U)
η_h = interpolate(η_cf,H)

f_scalar = panel_to_cartesian(f₀(ζ))
f_cf_scalar = panelwise_cellfield(f_scalar,Ω_panel,panel_ids)


tags = ["top_boundary", "bottom_boundary"]
Γ = BoundaryTriangulation(das,panel_model,tags=tags)
cell_geo_map = geo_map_func(get_panel_ids(Γ))
writevtk(Γ,dir*"/boundary",append=false,geo_map=cell_geo_map)

u_contra_skel = panelwise_cellfield(contra_v_3D(vX),Γ)
metric_skel = panelwise_cellfield(metric,Γ)
n_skel = CellField(VectorValue(1,0,0),Γ)

pts = get_cell_points(Γ)
# d = metric_skel⋅u_contra_skel
d = metric_skel
writevtk(Γ,dir*"/boundary",cellfields=["u"=>d],append=false,geo_map=cell_geo_map)



biformq(q,w) = ∫( (q⋅(metric_cf⋅w))*meas_cf )dΩ
liformq(w) = (
              ∫( (f_cf_scalar*(n⋅w))*(meas_cf/ff )  )dΩ
            + ∫( (metric_cf⋅u_contra_cf)⋅(curl(metric_cf⋅w) )  )dΩ
            # + ∫(  w⋅()  )dΛ
)
assem = SparseMatrixAssembler(H,R,das)
op = AffineFEOperator(biformq,liformq,H,R,assem)
qh = solve(ls,op)

qh_ambient = covarient_basis_cf ⋅ qh
η_ambient = covarient_basis_cf ⋅ η_h
e_η = l2((η_ambient - qh_ambient ),meas_cf,dΩ_error)

# latlon_cell_geo_map = geo_map_func(Ω_panel)
latlon_cell_geo_map = latlon_geo_map_func(Ω_panel)
cellfields = ["η"=>η_h, "ηambient"=>η_ambient,
              "η_h"=>qh ,"ηambient_h"=>qh_ambient,
              "eη_ambient"=>η_ambient-qh_ambient,"eη"=>η_h-qh,
              "uh"=>covarient_basis_cf⋅ u_contra_h]
 writevtk(Ω_panel,dir*"/latlon_solT_0.vtu", cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map)
