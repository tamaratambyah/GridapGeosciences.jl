module Distributed

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces
using LinearAlgebra
using FillArrays

using GridapDistributed
using P4est_wrapper
using GridapP4est
using PartitionedArrays

import GridapDistributed: DistributedDiscreteModel, GenericDistributedDiscreteModel
import GridapDistributed: BoundaryTriangulation
using GridapGeosciences.Geometry
import GridapGeosciences.Geometry: _CCAM_panel_wise_node_ids
import GridapGeosciences.Geometry: _CCAM_cube_nodes_3d
import GridapGeosciences.Geometry: setup_panel_cmaps
import GridapGeosciences.Geometry: panelwise_cellfield, geo_map_func, get_panel_ids

using GridapGeosciences.Fields
using GridapGeosciences.Visualisation

import Gridap.Visualization: writevtk, createvtk, write_vtk_file, create_vtk_file, create_pvtk_file

include("ParametricOctreeDistributedDiscreteModels.jl")
include("DistributedParametricDiscreteModel.jl")
include("panelwise_cellfield.jl")
include("panel_ids.jl")
include("Vtk.jl")
include("createpvd.jl")
include("helpers.jl")
include("BoundaryTriangulation.jl")

export ParametricOctreeDistributedDiscreteModel
export DistributedParametricDiscreteModel
export panelwise_cellfield, geo_map_func, get_panel_ids
export writevtk, createvtk, write_vtk_file, create_vtk_file, create_pvtk_file
export _make_pvd_distributed
export distributed_panel_ids
export DistributedAdaptivityGlue
export get_distributed_panel_model
export get_panel_ids, get_owned_panel_ids, get_skel_panel_ids
export BoundaryTriangulation

end
