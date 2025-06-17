using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays
import Gridap.TensorValues: meas
using StaticArrays

include("Geometry/coarse_cube_surface.jl")
include("Geometry/ManifoldGrid.jl")
include("Geometry/ManifoldDiscreteModel.jl")

include("Adaptivity/panel_ids_from_refinement.jl")
include("Adaptivity/Refinement.jl")

include("CellData/SurfaceMetric.jl")
include("CellData/SurfaceQuadrature.jl")

include("Operators/SurfaceOperators.jl")

include("Fields/Bump.jl")
include("Fields/PanelRotation.jl")
include("Fields/Gnomonic.jl")
include("Fields/Sigma.jl")
include("Fields/MetricMaps.jl")
include("Fields/FieldHelpers.jl")

include("FESpaces/FESpaceHelpers.jl")

include("Visualization/VisualizationData.jl")

include("ParametricAmbientSpaceMappers/VectorFields.jl")
include("ParametricAmbientSpaceMappers/ScalarFields.jl")


include("helpers.jl")


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

!isdir(plotsdir()) && mkdir(plotsdir())

function coarse_cube_model_3D(a::Real)

  cube_model_3D = _coarse_cube_model_3D(a)

  global A_bump, B_bump, b_bump = bump_matrics(a)
  global RADIUS = 1.0*sqrt(3.0)
  return cube_model_3D
end

function _coarse_cube_model_3D(a::Real)
  cube_grid_3D,topo_3D,face_labels_3D, = coarse_cube_surface_3D(a)
  cube_model_3D = UnstructuredDiscreteModel(cube_grid_3D,topo_3D,face_labels_3D)
  return cube_model_3D
end
