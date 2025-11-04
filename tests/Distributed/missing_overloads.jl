using Gridap
using GridapDistributed
using Gridap.FESpaces
using FillArrays


function Gridap.FESpaces._convert_to_collectable(object::Union{<:CellField,<:Function,<:Number},ntags)
  Gridap.FESpaces._convert_to_collectable(Fill(object,ntags),ntags)
end

function Gridap.FESpaces.TrialFESpace(f::GridapDistributed.DistributedSingleFieldFESpace,cf::GridapDistributed.DistributedCellField)
  println("my overload")
  spaces = map(f.spaces,cf.fields) do s, field
    TrialFESpace(s,field)
  end
  GridapDistributed.DistributedSingleFieldFESpace(spaces,f.gids,f.trian,f.vector_type,f.metadata)
end

function Gridap.FESpaces.interpolate_dirichlet!(
  u::GridapDistributed.DistributedCellField, free_values::AbstractVector,
  dirichlet_values::AbstractArray{<:AbstractVector},
  f::GridapDistributed.DistributedSingleFieldFESpace)
  println("interpolate dirichlt")
  map(local_views(u), f.spaces,local_views(free_values),dirichlet_values) do u,V,fvec,dvec
    interpolate_dirichlet!(u,fvec,dvec,V)
  end
  FEFunction(f,free_values,dirichlet_values)
end
