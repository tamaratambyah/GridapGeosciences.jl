module Adaptivity

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces
using LinearAlgebra
using FillArrays

using GridapGeosciences.Geometry
using GridapGeosciences.Fields
import GridapGeosciences.Geometry: CubedSphereParametricDiscreteModel, CubedSphere2DParametricDiscreteModel
import GridapGeosciences.Geometry: get_panel_ids, get_nodes_from_coords
import GridapGeosciences.Geometry: CubedSphereAmbientDiscreteModel
import GridapGeosciences.Fields: MyAffineField

include("Refinement.jl")

export refine
export is_child

end
