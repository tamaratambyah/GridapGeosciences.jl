module Distributed

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces
using LinearAlgebra
using FillArrays

using GridapDistributed
using GridapP4est
using PartitionedArrays

import GridapDistributed: DistributedDiscreteModel, GenericDistributedDiscreteModel

using GridapGeosciences.Geometry
import GridapGeosciences.Geometry: _CCAM_panel_wise_node_ids
import GridapGeosciences.Geometry: _CCAM_cube_nodes_3d

include("ParametricOctreeDistributedDiscreteModels.jl")

export ParametricOctreeDistributedDiscreteModel

end
