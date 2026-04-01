"""
    struct LinearCombinationDofVector{T<:Dof,V,F} <: AbstractVector{T}
      values :: V
      predofs:: F
    end

Type that implements a dof basis (a) as the linear combination of a dof pre-basis
(b). The dofs are first evaluated at the dof pre-basis (b) (field `predofs`) and the
predof values are next mapped to dof basis (a) applying a change of basis (field
`values`).

Fields:

- `values::AbstractMatrix{<:Number}` the matrix of the change from dof basis (b) to (a)
- `predofs::AbstractVector{T}` A type representing dof pre-basis (b), with `T<:Dof`
"""
struct LinearCombinationDofVector{T,V,F} <: AbstractVector{T}
  values::V
  predofs::F
  function LinearCombinationDofVector(
    values::AbstractMatrix{<:Number},
    predofs::AbstractVector{<:Dof}
  )
    @check size(values,1) == length(predofs) """\n
    Incompatible sizes for performing the linear combination

        linear_combination(values,predofs) = transpose(values)*predofs

    size(values,1) != length(predofs)
    """
    T = eltype(predofs)
    V = typeof(values)
    F = typeof(predofs)
    new{T,V,F}(values,predofs)
  end
end

Base.size(a::LinearCombinationDofVector) = (size(a.values,2),)
Base.IndexStyle(::LinearCombinationDofVector) = IndexLinear()
Base.getindex(::LinearCombinationDofVector{T},::Integer) where T = T()

function linear_combination(a::AbstractMatrix{<:Number}, b::AbstractVector{<:Dof})
  LinearCombinationDofVector(a,b)
end

function return_cache(b::LinearCombinationDofVector,field)
  k = LinearCombinationMap(:)
  cf = return_cache(b.predofs,field)
  fx = evaluate!(cf,b.predofs,field)
  ck = return_cache(k,fx,transpose(b.values))
  return cf, ck
end

function evaluate!(cache,b::LinearCombinationDofVector,field)
  cf, ck = cache
  k = LinearCombinationMap(:)
  fx = evaluate!(cf,b.predofs,field)
  return evaluate!(ck,k,fx,transpose(b.values))
end

# This constructor overload is NOT in Gridap 0.20.x. 
# We needed so that the method 
# Base.getindex(::LinearCombinationDofVector{T},::Integer) where T = T()
# above does not throw an error. In Gridap 0.20.x, there is NO constructor
# for PointValue without arguments. I think this should be reported as 
# an issue in Gridap.
function PointValue{Point{D,Float64}}() where D
   return PointValue(ntuple(i->0.0,D))
end