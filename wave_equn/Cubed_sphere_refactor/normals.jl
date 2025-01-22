function get_outward_normal_vector(model::CubedSphereDiscreteModel)
  cell_Jt, cell_x = get_tangent_space_basis(model)
  cell_normal = lazy_map( Operation(_unit_outward_normal), cell_Jt)
  Gridap.CellData.GenericCellField(cell_normal,Triangulation(model),ReferenceDomain())
end


"""
The columns of the Jacobian form a basis of the tangent space.
The tangent space is a vector space of dimension 2
"""
function get_tangent_space_basis(model::CubedSphereDiscreteModel)
  cell_reffe = get_reffes(get_grid(model))
  cell_dofs = lazy_map(get_dof_basis,cell_reffe)
  cell_x = lazy_map(get_nodes,cell_dofs)
  cell_map = get_cell_map(model)
  cell_Jt = lazy_map(Broadcasting(∇),cell_map)
  return cell_Jt, cell_x
end


"""
_unit_outward_normal -- takes the cross product of the basis vectors of the
tangent space, which is perpendicular to the tangent space. i.e. normal to the
surface

Note: the indexing of TensorValues is backwards to matrices. Index the "column"
first, and then the "row". A 2x3 TensorValue has the from
v = [1 2
     3 4
     5 6]
where
v[1,1] = 1;   v[2,1] = 2
v[1,2] = 3;   v[2,2] = 4
v[1,3] = 5;   v[2,3] = 6
"""
function _unit_outward_normal(v::Gridap.TensorValues.MultiValue{Tuple{2,3}})
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  n = VectorValue(n1,n2,n3)
  n/norm(n)
end
