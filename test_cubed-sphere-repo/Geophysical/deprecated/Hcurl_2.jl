using MPI
using PartitionedArrays

# add Gridap#cubed-sphere-moment-based-reffes

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
include(srcdir("Helpers/overloads.jl"))

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


#################### sphere
n_ref_lvls = 3

o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
        num_horizontal_uniform_refinements=n_ref_lvls,
        num_vertical_uniform_refinements=0)
panel_model = o3model.parametric_dmodel

das =  FullyAssembledRows()

panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(das,panel_model)
dΩ = Measure(Ω_panel,4)
dΩ_error = Measure(Ω_panel,8)
metric_cf = ParametricCellField(metric,Ω_panel,panel_ids)
covariant_basis_cf = ParametricCellField(covariant_basis,Ω_panel,panel_ids)

## finite element space with boundry conditions
tags = ["top_boundary", "bottom_boundary"]

p_fe = 1

# R = TestFESpace(panel_model, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
# H = TrialFESpace(R,VectorValue(1.0,0.0,0.0))
# f_cf = CellField(VectorValue(1.0,0.0,0.0),Ω_panel) ## normal vector in the chart

# f_vec₀(xyz) = VectorValue(xyz[1],xyz[2],0.0)

# function f_vec₀(xyz)
#     # f = η₀(0.0)(xyz)
#     # n = xyz[1]*normal_vec(xyz) # radial normal
#     # f*n
#     VectorValue(xyz[2]*xyz[1],xyz[2]*xyz[2],xyz[2]*xyz[3])
#     # VectorValue(xyz[1],xyz[2],xyz[3])
#     # n
# end

inv_jacobian(p) = x -> inv(forward_jacobian(p)(x))
contra_v_3D(vecX::Function,p::Int) = αβ -> inv_jacobian(p)(αβ) ⋅ vecX(p)(αβ)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)

# f_vec = panel_to_cartesian(f_vec₀)
# f_cf = ParametricCellField(contra_v_3D(f_vec),Ω_panel,panel_ids)

function fV(p)
  function f(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    VectorValue(xyz[1]*xyz[1],xyz[1]*xyz[2],xyz[1]*xyz[3])
  end
end


f_cf = ParametricCellField(contra_v_3D(fV),Ω_panel,panel_ids)

R = TestFESpace(panel_model, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl)
H = TrialFESpace(R)

f_h = interpolate(f_cf,H)
grad = gradient(f_cf)
gradh = gradient(f_h)

# check the XX component
scalar(xyz) = 2*xyz[1]
s_cf = ParametricCellField(panel_to_cartesian(scalar),Ω_panel,panel_ids)


latlon_cell_geo_map = latlon_geo_map_func(Ω_panel)
# latlon_cell_geo_map = geo_map_func(Ω_panel)
panel_cfs = [f_h, f_cf,  gradh, gradh, grad,  f_h-f_cf, gradh-grad, s_cf    ]
cellfields = map((x,y) -> x=>y, ["f_h", "f", "grad_h", "grad_h2", "grad",  "ef", "egrad" , "s"],panel_cfs)
writevtk(Ω_panel,dir*"/sol.vtu", cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map)
