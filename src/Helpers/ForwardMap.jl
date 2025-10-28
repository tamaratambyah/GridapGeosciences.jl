
using StaticArrays
struct ForwardMap{A}  <: Field
  p::A
end

function Gridap.Arrays.return_cache(f::ForwardMap,cellx::AbstractArray{<:VectorValue{2}})
  y = similar(cellx,VectorValue{3,Float64})
  return y
end


function Gridap.Arrays.evaluate!(cache,f::ForwardMap,cellx::AbstractArray{<:VectorValue{2}} )
  p = f.p
  y = cache

  map!(x-> forward_map_2D(p,x),
  y, cellx  )

  return y
end


function Gridap.Arrays.return_cache(f::ForwardMap,x::VectorValue{2})
  y = zero(VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::ForwardMap,x::VectorValue{2})
  p = f.p
  y = cache
  y = forward_map_2D(p,x)
  return y
end

"""
The jacobian is 3 x 2
  J = [dXda dXdb
       dYda dYdb
       dZda dZdb  ]
As a TensorValue data = (dXda,dYda,dZda,   dXdb, dYdb, dZdb)

Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)

The transpose is 2 x 3
  JT = [dXda dYda dZda
        dXdb dYdb dZdb  ]
As a TensorValue data = (dXda,dXdb,  dYda,dYdb,  dZda,dZdb)
"""
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMap},
  cellx::AbstractArray{<:VectorValue{2}})
  y = similar(cellx,TensorValue{2,3,Float64})
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:ForwardMap},
  cellx::AbstractArray{<:VectorValue{2}})
  cache, = c
  setsize!(cache,size(cellx))
  p = f.object.p

  y = cache
  # map!(x-> TensorValue{2,3}( dXda(x),dXdb(x),  dYda(x),dYdb(x),  dZda(x),dZdb(x) ),
      # y, cellx  )
  map!(x-> transpose(forward_jacobian_2D(p,x)),
      y, cellx  )
  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMap},x::VectorValue{2})
  zero(TensorValue{2,3,Float64})
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:ForwardMap},x::VectorValue{2})
  y = cache
  p = f.object.p

  y = transpose(forward_jacobian_2D(p,x))
  return y
end


################################################################################
########## 3D ########
################################################################################

function Gridap.Arrays.return_cache(f::ForwardMap,cellx::AbstractArray{<:VectorValue{3}})
  y = similar(cellx,VectorValue{3,Float64})
  return y
end


function Gridap.Arrays.evaluate!(cache,f::ForwardMap,cellx::AbstractArray{<:VectorValue{3}} )
  p = f.p
  y = cache
  map!(x-> forward_map_3D(p,x),
      y, cellx  )

  return y
end


function Gridap.Arrays.return_cache(f::ForwardMap,x::VectorValue{3})
  y = zero(VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::ForwardMap,x::VectorValue{3})
  p = f.p
  y = cache
  y =  forward_map_3D(p,x)
  return y
end

"""
The jacobian is 3 x 2
  J = [dXdg dXda dXdb
       dYdg dYda dYdb
       dZdg dZda dZdb  ]
As a TensorValue data = (dXdg,dYdg,dZdg,  dXda,dYda,dZda,   dXdb, dYdb, dZdb)

Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)

The transpose is 2 x 3
  JT = [dXdg dYdg dZdg
        dXda dYda dZda
        dXdb dYdb dZdb  ]
As a TensorValue data = (dXdg,dXda,dXdb,  dYdg,dYda,dYdb,  dZdg,dZda,dZdb)
"""
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMap},
  cellx::AbstractArray{<:VectorValue{3}})
  y = similar(cellx,TensorValue{3,3,Float64})
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:ForwardMap},
  cellx::AbstractArray{<:VectorValue{3}})
  cache, = c
  setsize!(cache,size(cellx))
  p = f.object.p

  y = cache
  # map!(x-> TensorValue{2,3}( dXda(x),dXdb(x),  dYda(x),dYdb(x),  dZda(x),dZdb(x) ),
      # y, cellx  )
  map!(x-> transpose(forward_jacobian_3D(p,x)),
      y, cellx  )
  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMap},x::VectorValue{3})
  zero(TensorValue{3,3,Float64})
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:ForwardMap},x::VectorValue{3})
  y = cache
  p = f.object.p
  y = transpose(forward_jacobian_3D(p,x))
  return y
end
