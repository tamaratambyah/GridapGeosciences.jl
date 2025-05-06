using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays
import Gridap.TensorValues: meas

a = 1.0

include("Geometry/coarse_cube_surface.jl")

include("Geometry/Bump.jl")
include("Geometry/PanelRotation.jl")

include("Geometry/ManifoldGrid.jl")
include("Geometry/ManifoldDiscreteModel.jl")

include("Adaptivity/panel_ids_from_refinement.jl")
include("Adaptivity/Refinement.jl")

include("helpers.jl")


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


cube_grid_3D,topo_3D,face_labels_3D, = coarse_cube_surface_3D(a)
cube_model_3D = UnstructuredDiscreteModel(cube_grid_3D,topo_3D,face_labels_3D)
