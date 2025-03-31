using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays


const a = π/4.0 ## to generate the cube of central angles
_a = 1.0 ## to compute the radius of the sphere
const r = _a*sqrt(3.0)



include(("../src/Maps/PanelRotation.jl"))
include("../src/Maps/Bump.jl")
include("../src/Maps/SphericalMaps.jl")

include("Geometry/cube_surface_1_cell_per_panel.jl")
include("Geometry/ManifoldGrid.jl")
include("Geometry/ManifoldDiscreteModel.jl")

include("Adaptivity/panel_ids_from_refinement.jl")
include("Adaptivity/Refinement.jl")

include("../src/helpers.jl")

# include("maps/metric_maps.jl")
# include("maps/cube_cell_map.jl")
# include("maps/sphere_cell_map.jl")

# include("surface_metric_and_op/metric_info.jl")
# include("surface_metric_and_op/cubedsphere_metric.jl")
# include("surface_metric_and_op/operators.jl")
# include("surface_metric_and_op/quadrature.jl")


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
