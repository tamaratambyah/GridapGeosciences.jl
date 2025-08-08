include("inverse_field.jl")

R1 = [1 0 0
      0 1 0
      0 0 1]
R2 = [0 0 1
      0 1 0
      1 0 0]
R3 = [0 1 0
      1 0 0
      0 0 1]
R4 = [-1 0 0
      0 -1 0
      0 0 1]
R5 = [0 0 -1
      0 1 0
      -1 0 0]
R6 = [0 -1 0
      -1 0 0
      0 0 1]
Rs = [R1,R2,R3,R4,R5,R6]

panel_grid = get_grid(panel_model)
panel_topo = get_grid_topology(panel_model)

cmaps = get_cell_map(panel_grid)

# forward_m = lazy_map(p->ForwardMapField(p),panel_ids)

forward_m = lazy_map(p->SwapField(TensorValue(Rs[p])) ∘ ForwardMapField(1),panel_ids)
Jt = lazy_map(Broadcasting(gradient),forward_m)


inv_f = lazy_map(p->InverseMapField(p),panel_ids)
ambient_maps = lazy_map(∘,forward_m,cmaps)

ref_points = get_cell_ref_coordinates(panel_grid)
ambient_cell_coords = lazy_map(evaluate,ambient_maps,ref_points)./1

function get_nodes_from_coords(grid::Grid{Dc,Dp},
  coords_array::AbstractArray{<:Vector{<:VectorValue{D,T}}}) where {Dc,Dp,D,T}

  cell_node_ids = get_cell_node_ids(grid)
  nodes = similar(coords_array, VectorValue{D,T}, num_nodes(grid))

  get_nodes_from_coords!(nodes,cell_node_ids,coords_array)

  return nodes
end

function get_nodes_from_coords!(nodes,cell_node_ids,
  coords_array::AbstractArray{<:Vector{<:VectorValue}})

  cache = array_cache(coords_array)

  for i in eachindex(coords_array)
    ids = cell_node_ids[i]
    nodes[ids] .= getindex!(cache, coords_array, i)
  end

end


ambient_nodes = get_nodes_from_coords(panel_grid,ambient_cell_coords)

ambient_grid = Geometry.UnstructuredGrid(ambient_nodes,get_cell_node_ids(panel_grid),get_reffes(panel_grid),get_cell_type(panel_grid),OrientationStyle(panel_grid),
                    nothing,ambient_maps)

ambient_topo = UnstructuredGridTopology(ambient_nodes,get_cell_node_ids(panel_grid),get_cell_type(panel_topo),get_polytopes(panel_topo),OrientationStyle(panel_topo))
ambient_labels = FaceLabeling(panel_topo)

ambient_model = UnstructuredDiscreteModel(ambient_grid,ambient_topo,ambient_labels)


ucf_mapped = lazy_map(Broadcasting(∘),get_data(ucf),inv_f)
ucf_ambient = CellData.GenericCellField(ucf_mapped,Triangulation(ambient_model),PhysicalDomain() )


_uh = change_domain(uh,ReferenceDomain(),PhysicalDomain())
uh_mapped = lazy_map(Broadcasting(∘),get_data(_uh),inv_f)
uh_ambient = CellData.GenericCellField(uh_mapped,Triangulation(ambient_model),PhysicalDomain() )

_uh_l2 = change_domain(uh_l2,ReferenceDomain(),PhysicalDomain())
uh_l2_mapped = lazy_map(Broadcasting(∘),get_data(_uh_l2),inv_f)
uh_l2_ambient = CellData.GenericCellField(uh_l2_mapped,Triangulation(ambient_model),PhysicalDomain() )


writevtk(Triangulation(ambient_model),dir*"/ambient_model",
cellfields=["u"=>ucf_ambient,"uh"=>uh_ambient, "uh_l2"=>uh_l2_ambient,"e"=>ucf_ambient-uh_ambient],
append=false)






mapping = map(p->SwapField(TensorValue(Rs[p])) ∘ ForwardMapField(1),panel_ids)
inv_mapping =  lazy_map(p->InverseMapField(p),panel_ids)

cell_field = map(p->GenericField(J(p)),panel_ids)
J_cf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

_cf = J_cf ⋅ gradient(ucf)
cf_mapped = lazy_map(Broadcasting(∘),get_data(_cf),inv_mapping)
cf_ambient = CellData.GenericCellField(cf_mapped,Triangulation(ambient_model),PhysicalDomain() ) # ambient cell field


cell_field = map(p->GenericField(Jcov(p)),panel_ids)
Jcov_cf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

s_grad = Jcov_cf ⋅ (inv_metric_cf ⋅ gradient(ucf))
s_grad_mapped = lazy_map(Broadcasting(∘),get_data(s_grad),inv_mapping)
s_grad_ambient = CellData.GenericCellField(s_grad_mapped,Triangulation(ambient_model),PhysicalDomain() )

writevtk(Triangulation(ambient_model),dir*"/ambient_model",
          cellfields=["gradu"=>cf_ambient,"u"=>ucf_ambient,"sgrad"=>s_grad_ambient],
          append=false)
