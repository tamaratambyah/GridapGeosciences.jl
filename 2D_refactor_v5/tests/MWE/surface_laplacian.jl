"""
Consider mapping a square x,y to another square x+2,y+2
  ϕ(x,y) = (x+2,y+2)
  This map is implemented as SquareShiftField

Consider mapping a square x,y to a plane via a quadratic mapping:
  φ(x,y) = (x,y,x^2 + y^2)
  This map is implemented as QuadPlaneField
This mapping has a non-constant metric.

Consider the composition: φ ∘ ϕ
This take a point in 2D parametric space -> 3D ambient space

Consider a function defined in the ambient coordinate system XYZ.
Want to map this function back to the square, and compute the surface laplacian
"""

using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields, Gridap.TensorValues
using StaticArrays


"""
SquareShiftField: 2D->2D
ϕ(x,y) = (x+2,y+2)
"""
struct SquareShiftField <: Gridap.Fields.Field
end

function Gridap.Arrays.return_cache(k::SquareShiftField,cellx::AbstractArray{<:VectorValue{2}})
  y = similar(cellx,VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SquareShiftField,cellx::AbstractArray{<:VectorValue{2}})
  # setsize!(cache,size(cellx))
  y = cache
  map!(x -> VectorValue(x[1]+2, x[2]+2),
                        y, cellx)
  return y

end

function Gridap.Arrays.return_cache(k::SquareShiftField,x::VectorValue{2})
  y = zero(VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SquareShiftField,x::VectorValue{2})
  y = cache
  y = VectorValue(x[1]+2, x[2]+2)
  return y
end

function Gridap.Arrays.return_cache(k::SquareShiftField,x::SVector{2})
  y = zero(VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SquareShiftField,x::SVector{2})
  y = cache
  y = VectorValue(x[1]+2, x[2]+2)
  return y
end

"""
gradient of the quadratic mapping: (x,y) → (x+2,y+2) is:
  J = [ 1 0
        0 1 ]
"""

function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:SquareShiftField},cellx::AbstractArray{<:VectorValue{2}})
  _T = typeof(TensorValue{2,2,Float64})
  y = similar(cellx,_T,size(cellx))
  CachedArray(y)
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:SquareShiftField},cellx::AbstractArray{<:VectorValue{2}})
  cache,  = c
  setsize!(cache,size(cellx))
  y = cache.array

  map!(x -> TensorValue{2,2}( 1.0,0.0,  0.0,1.0 ),
                    y, cellx)

  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:SquareShiftField},x::VectorValue{2})
  zero(TensorValue{2,2,Float64})
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:SquareShiftField},x::VectorValue{2})
  y = cache
  y = TensorValue{2,3}( 1.0,0.0,  0.0,1.0 )
  return y
end



"""
QuadPlaneField: 2D-> 3D
φ(x,y) = (x,y,x^2 + y^2)
"""
struct QuadPlaneField <: Gridap.Fields.Field
end

function Gridap.Arrays.return_cache(k::QuadPlaneField,cellx::AbstractArray{<:VectorValue{2}})
  y = similar(cellx,VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::QuadPlaneField,cellx::AbstractArray{<:VectorValue{2}})
  # setsize!(cache,size(cellx))
  y = cache
  map!(x -> VectorValue(x[1], x[2], x[1]^2 + x[2]^2),
                        y, cellx)
  return y

end

function Gridap.Arrays.return_cache(k::QuadPlaneField,x::VectorValue{2})
  y = zero(VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::QuadPlaneField,x::VectorValue{2})
  y = cache
  y = VectorValue(x[1], x[2], x[1]^2 + x[2]^2)
  return y
end

function Gridap.Arrays.return_cache(k::QuadPlaneField,x::SVector{2})
  y = zero(VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::QuadPlaneField,x::SVector{2})
  y = cache
  y = VectorValue(x[1], x[2], x[1]^2 + x[2]^2)
  return y
end



"""
gradient of the quadratic mapping: (x,y) → XYZ is:
  J = [ 1 0
        0 1
        2x 2y ]
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
# Note  TensorValue{2,3}(0,4,1,0,5,1) == [0 1 5
                                          4 0 1]
The transpose of the jacobian is:
  JT = [ 1 0 2x
         0 1 2y]
To write as a TensorValue:
  TensorValue{2,3}(1,0,0,1,2x,2y)
"""

function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:QuadPlaneField},cellx::AbstractArray{<:VectorValue{2}})
  _T = typeof(TensorValue{2,3,Float64})
  y = similar(cellx,_T,size(cellx))
  CachedArray(y)
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:QuadPlaneField},cellx::AbstractArray{<:VectorValue{2}})
  cache,  = c
  setsize!(cache,size(cellx))
  y = cache.array

  map!(x -> TensorValue{2,3}( 1.0,0.0,  0.0,1.0,  2*x[1],2*x[2]  ),
                    y, cellx)

  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:QuadPlaneField},x::VectorValue{2})
  zero(TensorValue{2,3,Float64})
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:QuadPlaneField},x::VectorValue{2})
  y = cache
  y = TensorValue{2,3}( 1.0,0.0,  0.0,1.0,  2*x[1],2*x[2] )
  return y
end


##### Apply mappings

model_parametric = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(2,2)))
Ω_parametric = Triangulation(model_parametric)
dΩ_parametric = Measure(Ω_parametric,2)
pt = Point(0.0,0.0)


### quadratic mapping, φ(x,y) = (x,y,x^2 + y^2)
metric_func(x) = TensorValue{2,2}(1+4*x[1]^2,4*x[1]*x[2],4*x[1]*x[2],1+4*x[2]^2 )
sq_meas_func(x) = sqrt(meas(metric_func(x)))
inv_metric_func(x) = inv(metric_func(x))

m_cf = CellField(metric_func,Ω_parametric) # 2x2 tensor
sq_meas = CellField(sq_meas_func,Ω_parametric ) # scalar
inv_metric = CellField(inv_metric_func,Ω_parametric) # 2x2 tensor

function fambient_scalar(X)
  X[1]*X[2]*X[3]
end

################################################################################
#### Consider mapping via separate evaluation of the fields (i.e. no composition)
#### This implementation works
################################################################################
function fparametric_scalar(xy)
  _xy = SquareShiftField()(xy)
  XYZ = QuadPlaneField()(_xy)
  fambient_scalar(XYZ)
end

cellf = map(p->GenericField(fparametric_scalar), collect(1:num_cells(model_parametric)))
cf_parametric = CellData.GenericCellField(cellf,Ω_parametric,PhysicalDomain())

### laplacian is the divergence(gradient)
surfgrad = inv_metric ⋅ gradient(cf_parametric)
f = sq_meas*surfgrad
surflap = 1/sq_meas * divergence(f)

surflap(pt)

################################################################################
#### Consider mapping via composition
#### This implementation fails
################################################################################
function _fparametric_scalar(xy)
  cmap = QuadPlaneField() ∘  SquareShiftField()
  XYZ = evaluate(cmap,xy)

  fambient_scalar(XYZ)
end

cellf = map(p->GenericField(_fparametric_scalar), collect(1:num_cells(model_parametric)))
cf_parametric = CellData.GenericCellField(cellf,Ω_parametric,PhysicalDomain())

### laplacian is the divergence(gradient)
_surfgrad = inv_metric ⋅ gradient(cf_parametric) ## FAILS!
_f = sq_meas*_surfgrad
_surflap = 1/sq_meas * divergence(f_)

_surflap(pt)
