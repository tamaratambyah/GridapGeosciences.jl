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

dir = datadir("Hcurl")
!isdir(dir) && mkdir(dir)

include(srcdir("Helpers/overloads.jl"))
include("../../Geophysical/CurlConformingFESpacesFixes.jl")

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


#################### sphere
o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
        num_horizontal_uniform_refinements=0,
        num_vertical_uniform_refinements=0)
panel_model = o3model.parametric_dmodel




tags = ["top_boundary", "bottom_boundary"]

panel_ids = get_panel_ids(panel_model)
Ω = Triangulation(panel_model)
dΩ = Measure(Ω,6)
dΩ_error = Measure(Ω,8)
