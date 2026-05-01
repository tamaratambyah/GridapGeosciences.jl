module ConvergenceTools


using Test
using Gridap
using Gridap.Geometry, Gridap.Adaptivity
using GridapDistributed
using GridapP4est
using PartitionedArrays

using GridapGeosciences.Geometry
using GridapGeosciences.Distributed

include("convergence_tools.jl")

export p_convergence_auto_test, h_convergence_auto_test
export get_distributed_refined_models, get_octree_refined_models, get_3D_octree_refined_models
export nref, nc, nc_horizontal, dx, dx_horizontal
export convergence_rate


end
