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
include(srcdir("Helpers/overloads.jl"))
include("CurlConformingFESpacesFixes.jl")

## pullback 3D vector to 3D chart
inv_jacobian(p) = x -> inv(forward_jacobian_3D(p)(x))
contra_v_3D(vecX::Function,p::Int) = αβ -> inv_jacobian(p)(αβ) ⋅ vecX(p)(αβ)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


#################### sphere
n_ref_lvls = 2

o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
        num_horizontal_uniform_refinements=n_ref_lvls,
        num_vertical_uniform_refinements=0)
panel_model = o3model.parametric_dmodel


panel_ids = get_panel_ids(panel_model)
Ω = Triangulation(panel_model)
dΩ = Measure(Ω,4)

## normal vector in the chart
# f_cf = CellField(VectorValue(1,0.0,0.0),Ω)
function fV(p)
  function f(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    VectorValue(xyz[1],xyz[2],xyz[3])
  end
end
f_cf = panelwise_cellfield(contra_v_3D(fV),Ω,panel_ids)

## finite element space with boundry conditions
tags = ["top_boundary", "bottom_boundary"]
p_fe = 1

R = TestFESpace(panel_model, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
H = TrialFESpace(R,f_cf)

f_h = interpolate(f_cf,H)
eh = f_cf-f_h


print("error", sum(∫( eh⋅eh )*dΩ)); print("\n");
@assert sum(∫( eh⋅eh )*dΩ)/sum(∫( f_cf⋅f_cf )*dΩ) < 1.e-12

grad = gradient(f_cf)
gradh = gradient(f_h)

eh_grad = grad-gradh
@assert sum(∫( eh_grad⊙eh_grad )*dΩ) < 1.e-12

latlon_geo_map = latlon_geo_map_func(Ω)
panel_cfs = [f_h, f_cf, eh, gradh, grad, eh_grad   ]
cellfields = map((x,y) -> x=>y, ["f_h", "f", "eh", "grad_h", "grad", "eh_grad" ],panel_cfs)
writevtk(Ω,dir*"/sol.vtu", cellfields=cellfields,append=false,geo_map=latlon_geo_map)


################################################################################
#### Consider a more complicated vector field
################################################################################
function fV(p)
  function f(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    VectorValue(-xyz[2],xyz[1],0)
  end
end

f_cf = panelwise_cellfield(contra_v_3D(fV),Ω,panel_ids)
H = TrialFESpace(R,f_cf)

f_h = interpolate(f_cf,H)
grad = gradient(f_cf)
gradh = gradient(f_h)

dΩ_error = Measure(Ω,6)

eh = f_cf-f_h
sum(∫( eh⋅eh )*dΩ_error)
sum(∫( eh⋅eh )*dΩ_error)/sum(∫( f_cf⋅f_cf )*dΩ_error)

eh_grad = grad-gradh
sum(∫( eh⊙eh )*dΩ_error)

latlon_geo_map = latlon_geo_map_func(Ω)
panel_cfs = [f_h, f_cf, eh, gradh, grad, eh_grad   ]
cellfields = map((x,y) -> x=>y, ["f_h", "f", "ef", "grad_h", "grad", "egrad" ],panel_cfs)
writevtk(Ω,dir*"/sol.vtu", cellfields=cellfields,append=false,geo_map=latlon_geo_map)
