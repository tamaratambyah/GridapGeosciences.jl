using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using LinearAlgebra
using FillArrays

model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(3,3)))

topo = get_grid_topology(model)

Table(get_cell_permutations(topo))
