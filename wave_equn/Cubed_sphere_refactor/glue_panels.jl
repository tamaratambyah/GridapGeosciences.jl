using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Test
using LinearAlgebra
using FillArrays
include("individual_panels.jl")


ncells_per_panel = 4
npanels = 6
ref_panel = CartesianDiscreteModel((-1,1,-1,1),(ncells_per_panel,ncells_per_panel))

cell_coordinates = get_cell_coordinates(ref_panel)
cell_node_ids = get_cell_node_ids(ref_panel)
node_coordinates = get_node_coordinates(ref_panel)

panel1 = get_panel(1,node_coordinates,cell_node_ids)
panel2 = get_panel(2,node_coordinates,cell_node_ids)



function generate_ptr(n)
  nvertices = 4
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end

ptr  = [ 1, 5, 9, 13, 17, 21, 25 ]
ptr = generate_ptr(npanels)
data = [ 2,4,6,8, 4,3,8,7, 3,1,7,5, 1,2,5,6, 6,8,5,7, 1,3,2,4  ]
cube_vertex_ids = Gridap.Arrays.Table(data,ptr)
