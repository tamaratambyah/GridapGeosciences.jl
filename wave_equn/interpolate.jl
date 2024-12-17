"""
function interpolate(object, fs::SingleFieldFESpace)
  free_values = zero_free_values(fs)
  interpolate!(object, free_values,fs)
end

function interpolate!(object, free_values,fs::SingleFieldFESpace)
  cell_vals = _cell_vals(fs,object)
  gather_free_values!(free_values,fs,cell_vals)
  FEFunction(fs,free_values)
end

function _cell_vals(fs::SingleFieldFESpace,object)
  s = get_fe_dof_basis(fs)
  trian = get_triangulation(s)
  f = CellField(object,trian,DomainStyle(s))
  cell_vals = s(f)
end

function gather_free_values!(free_values,f::SingleFieldFESpace,cell_vals)
    dirichlet_values = zero_dirichlet_values(f)
    gather_free_and_dirichlet_values!(free_values,dirichlet_values,f,cell_vals)
    free_values
end

function gather_free_and_dirichlet_values!(free_vals,dirichlet_vals,f::UnconstrainedFESpace,cell_vals)

  cell_dofs = get_cell_dof_ids(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  cells = 1:length(cell_vals)

  _free_and_dirichlet_values_fill!(
    free_vals,
    dirichlet_vals,
    cache_vals,
    cache_dofs,
    cell_vals,
    cell_dofs,
    cells)

  (free_vals,dirichlet_vals)
end

function  _free_and_dirichlet_values_fill!(
  free_vals,
  dirichlet_vals,
  cache_vals,
  cache_dofs,
  cell_vals,
  cell_dofs,
  cells)

  for cell in cells
    vals = getindex!(cache_vals,cell_vals,cell)
    dofs = getindex!(cache_dofs,cell_dofs,cell)
    for (i,dof) in enumerate(dofs)
      val = vals[i]
      if dof > 0
        free_vals[dof] = val
      elseif dof < 0
        dirichlet_vals[-dof] = val
      else
        @unreachable "dof ids either positive or negative, not zero"
      end
    end
  end

end

function FEFunction(fe::SingleFieldFESpace, free_values)
  diri_values = get_dirichlet_dof_values(fe)
  FEFunction(fe,free_values,diri_values)
end


function FEFunction(
  fs::SingleFieldFESpace, free_values::AbstractVector, dirichlet_values::AbstractVector)
  cell_vals = scatter_free_and_dirichlet_values(fs,free_values,dirichlet_values)
  cell_field = CellField(fs,cell_vals)
  SingleFieldFEFunction(cell_field,cell_vals,free_values,dirichlet_values,fs)
end


function scatter_free_and_dirichlet_values(f::UnconstrainedFESpace,free_values,dirichlet_values)
  @check eltype(free_values) == eltype(dirichlet_values) "
  cell_dof_ids = get_cell_dof_ids(f)
  lazy_map(Broadcasting(PosNegReindex(free_values,dirichlet_values)),cell_dof_ids)
end

"""
