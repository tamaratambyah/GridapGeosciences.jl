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
include("Williamson2Test.jl")

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

n_ref_lvls = 3
ζ = 0.0

function η_vec₀(ζ)
  function _η₀(xyz)
    η = η₀(ζ)(xyz)
    n = normal_vec(xyz)
    η*n
  end
end

function f_vec₀(ζ)
  function _f₀(xyz)
    f = f₀(ζ)(xyz)
    n = normal_vec(xyz)
    f*n
  end
end


vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
f = panel_to_cartesian(f₀(ζ))
f_vec = panel_to_cartesian(f_vec₀(ζ))
η = panel_to_cartesian(η₀(ζ))
η_vec = panel_to_cartesian(η_vec₀(ζ))


o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
        num_horizontal_uniform_refinements=n_ref_lvls,
        num_vertical_uniform_refinements=0)
panel_model = o3model.parametric_dmodel

ls = LUSolver()
p_fe = 1

das =  FullyAssembledRows()

## finite element solver
panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(das,panel_model)
dΩ = Measure(Ω_panel,4*(p_fe+1))
dΩ_error = Measure(Ω_panel,8*(p_fe+1))

# Λ = SkeletonTriangulation(das,panel_model)
# dΛ = Measure(Λ,4*(p_fe+1))

tags = ["top_boundary", "bottom_boundary"]

V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

R = TestFESpace(panel_model, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
H = TrialFESpace(R,VectorValue(0.0,0.0,0.0))

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


# initial conditions
u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
u_contra_h = interpolate(u_contra_cf,U)
cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
f_cf = panelwise_cellfield(contra_v(f_vec),Ω_panel,panel_ids)
f_ambient = cor_cf*n_ambient

η_cf = panelwise_cellfield(η,Ω_panel,panel_ids)
η_ambient = η_cf*n_ambient
### put vector into FE space
_η_cf = panelwise_cellfield(contra_v(η_vec),Ω_panel,panel_ids)
η_h = interpolate(_η_cf,H)
_η_ambient = covarient_basis_cf ⋅ η_h
l2(η_ambient-_η_ambient,meas_cf,dΩ_error)


# plot
# latlon_cell_geo_map = geo_map_func(Ω_panel)
latlon_cell_geo_map = latlon_geo_map_func(Ω_panel)
panel_cfs = [f_ambient,n_ambient,η_ambient, _η_ambient, η_ambient-_η_ambient  ]
cellfields = map((x,y) -> x=>y, [ "f", "n", "eta1", "eta2","eeta"],panel_cfs)
writevtk(Ω_panel,dir*"/latlon_solT_0.vtu", cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map)


Aperp = [0 0 0
        0 0 -1
        0 1 0]
Rperp = TensorValue(Aperp)
Rperp_cf = CellField(Rperp,Ω_panel)

biformq(q,w) = ∫( (q⋅(metric_cf⋅w))*meas_cf )dΩ
liformq(w) = (
              ∫( (f_cf⋅(metric_cf⋅w))*(meas_cf )  )dΩ
              # ∫( cor_cf*(w⋅n)*(meas_cf/ff )  )dΩ
            + ∫( (metric_cf⋅u_contra_h)⋅(curl(metric_cf⋅w) )  )dΩ
)
assem = SparseMatrixAssembler(H,R,das)
op = AffineFEOperator(biformq,liformq,H,R,assem)
qh = solve(ls,op)

qh_ambient = covarient_basis_cf ⋅ qh
e_η = l2((_η_ambient - qh_ambient ),meas_cf,dΩ_error)


# latlon_cell_geo_map = geo_map_func(Ω_panel)
latlon_cell_geo_map = latlon_geo_map_func(Ω_panel)
panel_cfs = [ η_h, _η_ambient, qh ,qh_ambient, _η_ambient-qh_ambient ]
cellfields = map((x,y) -> x=>y, ["η", "ηambient", "η_h", "ηambient_h", "eη"  ],panel_cfs)
writevtk(Ω_panel,dir*"/latlon_solT_0.vtu", cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map)
