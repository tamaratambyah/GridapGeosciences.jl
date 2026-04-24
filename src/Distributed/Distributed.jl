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

import GridapDistributed: DistributedCellField, DistributedTriangulation
import GridapDistributed: DistributedFaceLabeling

import GridapDistributed: DistributedDiscreteModel, GenericDistributedDiscreteModel
import GridapDistributed: BoundaryTriangulation
using GridapGeosciences.Geometry
import GridapGeosciences.Geometry: _CCAM_panel_wise_node_ids
import GridapGeosciences.Geometry: _CCAM_cube_nodes_3d
import GridapGeosciences.Geometry: setup_panel_cmaps
import GridapGeosciences.Geometry: panelwise_cellfield, geo_map_func, latlon_geo_map_func
import GridapGeosciences.Geometry: get_panel_ids, get_forward_map_generator, get_radius, get_thickness
import GridapGeosciences.Geometry: pullback_area_form, pushforward_normal
import GridapGeosciences.Fields: ForwardMap
using GridapGeosciences.Fields
using GridapGeosciences.Visualisation
import GridapGeosciences.Fields: Cartesian2SphericalMap

import Gridap.Visualization: writevtk, createvtk, write_vtk_file, create_vtk_file, create_pvtk_file

include("ParametricOctreeDistributedDiscreteModels.jl")
include("Parametric3DOctreeDistributedDiscreteModels.jl")
include("DistributedParametricDiscreteModel.jl")
include("panelwise_cellfield.jl")
include("panel_ids.jl")
include("Vtk.jl")
include("createpvd.jl")
include("helpers.jl")
include("Triangulations.jl")

import Gridap.FESpaces: FESpace
import GridapDistributed: generate_gids, _find_vector_type, _add_distributed_constraint, DistributedSingleFieldFESpace, _setup_dmodel_and_dtrian
import GridapGeosciences.FESpaces: _generate_face_to_master_cell_id, _generate_change_of_basis_matrices, _get_value_type, ParamTrianType
include("GradConformingFESpaces.jl")

export ParametricOctreeDistributedDiscreteModel
export Parametric3DOctreeDistributedDiscreteModel
export DistributedParametricDiscreteModel
export panelwise_cellfield, geo_map_func, get_panel_ids, latlon_geo_map_func
export get_forward_map_generator, get_radius, get_thickness
export writevtk, createvtk, write_vtk_file, create_vtk_file, create_pvtk_file
export _make_pvd_distributed
export distributed_panel_ids
export DistributedAdaptivityGlue
export get_distributed_panel_model
export get_panel_ids, get_owned_panel_ids, get_skel_panel_ids
# export BoundaryTriangulation
export pullback_area_form
export pushforward_normal

end
