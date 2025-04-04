using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays
import Gridap.TensorValues: meas

const a = π/4.0 ## to generate the cube of central angles
_a = 1.0 ## to compute the radius of the sphere
const r = _a*sqrt(3.0)



include("Fields/PanelRotation.jl")
include("Fields/Bump.jl")
include("Fields/Gnomonic.jl")
include("Fields/Sigma.jl")
include("Fields/CentralAngle.jl")
include("Fields/MetricMaps.jl")

include("Geometry/cube_surface_1_cell_per_panel.jl")
include("Geometry/ManifoldGrid.jl")
include("Geometry/ManifoldDiscreteModel.jl")

include("Adaptivity/panel_ids_from_refinement.jl")
include("Adaptivity/Refinement.jl")

include("CellData/SurfaceMetric.jl")
include("CellData/SurfaceOperators.jl")
include("CellData/SurfaceQuadrature.jl")

include("helpers.jl")

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


cube_model_3D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel(a)...)

##### make some refined models
# model = cube_model_3D
# ref_model = Gridap.Adaptivity.refine(model)
# ref_ref_model = Gridap.Adaptivity.refine(ref_model)
# ref_ref_ref_model = Gridap.Adaptivity.refine(ref_ref_model)
# ref_ref_ref_ref_model = Gridap.Adaptivity.refine(ref_ref_ref_model)
# ref_ref_ref_ref_ref_model = Gridap.Adaptivity.refine(ref_ref_ref_ref_model)
