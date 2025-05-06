using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays
import Gridap.TensorValues: meas

include("Geometry/coarse_cube_surface.jl")
include("Geometry/ManifoldGrid.jl")
include("Geometry/ManifoldDiscreteModel.jl")

include("Adaptivity/panel_ids_from_refinement.jl")
include("Adaptivity/Refinement.jl")

include("CellData/SurfaceMetric.jl")
include("CellData/SurfaceOperators.jl")
include("CellData/SurfaceQuadrature.jl")

include("Fields/Bump.jl")
include("Fields/PanelRotation.jl")
include("Fields/MetricMaps.jl")

include("helpers.jl")


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

function coarse_cube_model_3D(a::Real)
  cube_grid_3D,topo_3D,face_labels_3D, = coarse_cube_surface_3D(a)
  cube_model_3D = UnstructuredDiscreteModel(cube_grid_3D,topo_3D,face_labels_3D)

  global A_bump, B_bump, b_bump = bump_matrics(a)

  return cube_model_3D
end
