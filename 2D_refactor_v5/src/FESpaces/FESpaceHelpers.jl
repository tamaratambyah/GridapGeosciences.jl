"""
_cell_vals
evaluate cell-wise arry of generic fields on the dof basis
"""
function Gridap.FESpaces._cell_vals(fs::SingleFieldFESpace,object::AbstractArray{<:GenericField})
  println("my cell vals")
  s = get_fe_dof_basis(fs)
  trian = get_triangulation(s)
  f = CellData.GenericCellField(object,trian,PhysicalDomain())
  s(f)
end

"""
_cell_vals
evaluate cell field on the dof basis
"""
function Gridap.FESpaces._cell_vals(fs::SingleFieldFESpace,cf_ambient::CellField)
  println("my cell vals")
  s = get_fe_dof_basis(fs)
  s(cf_ambient)
end
