using Gridap
include("../../src/initialise.jl")

"""
Test rotation of panels in 2D and 3D for the coarsest model with 1 cell per panel
"""

panel_ids = 1:6

################################################################################
#### Test 2D rotations
################################################################################
data_2D = [ 1,2,3,4, 1,2,3,4, 2,4,1,3, 4,3,2,1, 4,3,2,1, 2,4,1,3 ] # orientation of panels in 2D
ptr = generate_ptr(6)
cell_node_ids_2D = Table(data_2D,ptr)

nodes_2D = [
    Point(-1.0, -1.0)  # node 1
    Point(1.0, -1.0)   # node 2
    Point(-1.0, 1.0)   # node 3
    Point(1.0, 1.0)    # node 4
  ]
cell_nodes_2D = lazy_map(Broadcasting(Reindex(nodes_2D)),cell_node_ids_2D)

# map panel 1 -> panel p in 2D
cell_nodes_panel1 = fill(nodes_2D,6)
mapped_cell_nodes = lazy_map(R1pPanelMap2D(),cell_nodes_panel1,panel_ids)

# check mapped nodes are equiv to reindexed nodes
@test cell_nodes_2D == mapped_cell_nodes

# map panel p -> panel 1 in 3D
_cell_nodes_panel1 = lazy_map(Rp1PanelMap2D(),cell_nodes_2D,panel_ids)
for i = 1:6
  @test _cell_nodes_panel1[i] == cell_nodes_2D[1] # check each panel maps back to panel 1
end



################################################################################
#### Test 3D rotations
################################################################################
cell_node_ids_3D = get_cell_node_ids(get_grid(cube_model_2D))
nodes_3D = [
    Point(1.0, -1.0, -1.0)  # node 1
    Point(1.0, 1.0, -1.0)   # node 2
    Point(1.0, -1.0, 1.0)   # node 3
    Point(1.0, 1.0, 1.0)    # node 4
    Point(-1.0, -1.0, 1.0)  # node 5
    Point(-1.0, 1.0, 1.0)   # node 6
    Point(-1.0, 1.0, -1.0)  # node 7
    Point(-1.0, -1.0, -1.0) # node 8
  ]
cell_nodes_3D = lazy_map(Broadcasting(Reindex(nodes_3D)),cell_node_ids_3D)

# map panel 1 -> panel p in 3D
cell_nodes_panel1 = fill(cell_nodes_3D[1],6)
mapped_cell_nodes = lazy_map(R1pPanelMap3D(),cell_nodes_panel1,panel_ids)

# check mapped nodes are equiv to reindexed nodes
@test cell_nodes_3D == mapped_cell_nodes

# map panel p -> panel 1 in 3D
_cell_nodes_panel1 = lazy_map(Rp1PanelMap3D(),cell_nodes_3D,panel_ids)
for i = 1:6
  @test _cell_nodes_panel1[i] == cell_nodes_3D[1] # check each panel maps back to panel 1
end
