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

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


#################### sphere
n_ref_lvls = 2

o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
        num_horizontal_uniform_refinements=n_ref_lvls,
        num_vertical_uniform_refinements=0)
panel_model = o3model.parametric_dmodel

das =  FullyAssembledRows()

panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(das,panel_model)
dΩ = Measure(Ω_panel,4)
dΩ_error = Measure(Ω_panel,8)

## finite element space with boundry conditions
tags = ["top_boundary", "bottom_boundary"]

p_fe = 1

R = TestFESpace(panel_model, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
# R = TestFESpace(panel_model, ReferenceFE(lagrangian,VectorValue{3,Float64},p_fe);conformity=:H1,dirichlet_tags=tags)
H = TrialFESpace(R,VectorValue(0.0,0.0,0.0))

## normal vector in the chart
f_cf = CellField(VectorValue(1,0.0,0.0),Ω_panel)
f_h = interpolate(f_cf,H)

grad = gradient(f_cf)
gradh = gradient(f_h)

latlon_cell_geo_map = latlon_geo_map_func(Ω_panel)
panel_cfs = [f_h, f_cf,  gradh, grad,    ]
cellfields = map((x,y) -> x=>y, ["f_h", "f", "grad_h", "grad" ],panel_cfs)
writevtk(Ω_panel,dir*"/sol.vtu", cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map)
