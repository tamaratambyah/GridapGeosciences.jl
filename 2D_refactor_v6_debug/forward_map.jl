function forward_map(αβ,p)
  α,β = αβ


  rho = sqrt(1 + (tan(α))^2 + (tan(β))^2 )

  if p == 1
    X = 1/rho
    Y = 1/rho * tan(α)
    Z = 1/rho * tan(β)
  elseif p == 2
    X = -1/rho * tan(β)
    Y = 1/rho * tan(α)
    Z = 1/rho
  elseif p == 3
    X = -1/rho * tan(α)
    Y = 1/rho
    Z = 1/rho * tan(β)
  elseif p == 4
    X = -1/rho
    Y = 1/rho * tan(β)
    Z = 1/rho * tan(α)
  elseif p == 5
    X = -1/rho * tan(α)
    Y = 1/rho * tan(β)
    Z = -1/rho
  elseif p == 6
    X = -1/rho * tan(β)
    Y = -1/rho
    Z = 1/rho * tan(α)
  end


  RADIUS*Point(X,Y,Z)
end

function forward_jacobian(αβ,p)

  α,β = αβ
  rho = sqrt(1 + (tan(α))^2 + (tan(β))^2 )
  drho_da = - tan(α)*(sec(α))^2 / ( rho^3  )
  drho_db = - tan(β)*(sec(β))^2 / ( rho^3 )

  if p == 1
    dXda = drho_da
    dXdb = drho_db
    dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
    dYdb = drho_db*tan(α)
    dZda = drho_da*tan(β)
    dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2
  elseif p == 2
    dXda = -( drho_da*tan(β) )
    dXdb = -( drho_db*tan(β) + 1/rho*(sec(β))^2 )
    dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
    dYdb = drho_db*tan(α)
    dZda = drho_da
    dZdb = drho_db
  elseif p == 3
    dXda = -( drho_da*tan(α) + 1/rho*(sec(α))^2  )
    dXdb = -( drho_db*tan(α) )
    dYda = drho_da
    dYdb = drho_db
    dZda = drho_da*tan(β)
    dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2
  elseif p == 4
    dXda = -drho_da
    dXdb = -drho_db
    dYda = -(-drho_da*tan(β) )
    dYdb = -(-drho_db*tan(β) + -1/rho*(sec(β))^2 )
    dZda = -(-drho_da*tan(α) + -1/rho*(sec(α))^2  )
    dZdb = -(-drho_db*tan(α) )
  elseif p == 5
    dXda = -drho_da*tan(α) + -1/rho*(sec(α))^2
    dXdb = -drho_db*tan(α)
    dYda = -(-drho_da*tan(β))
    dYdb = -(-drho_db*tan(β) + -1/rho*(sec(β))^2 )
    dZda = -drho_da
    dZdb = -drho_db
  elseif p == 6
    dYda = -drho_da
    dYdb = -drho_db
    dZda = -(-drho_da*tan(α) + -1/rho*(sec(α))^2  )
    dZdb = -(-drho_db*tan(α) )
    dXda = -drho_da*tan(β)
    dXdb = -drho_da*tan(β) + -1/rho*(sec(β))^2
  end

  ## J = [dXda dXdb
  ##      dYda dYdb
  ##      dZda dZdb  ]
  ## As a TensorValue data = (dXda,dYda,dZda, dXdb, dYdb, dZdb)

  RADIUS*TensorValue{3,2}(dXda,dYda,dZda, dXdb,dYdb,dZdb)

end

function forward_jacobian(p)
  function _forward_jacobian(αβ)
    forward_jacobian(αβ,p)
  end
end

function metric(p)
  function _metric(αβ)
    transpose(forward_jacobian(αβ,p)) ⋅ forward_jacobian(αβ,p)
  end
end

function measure(p)
  function _measure(αβ)
    g = transpose(forward_jacobian(αβ,p)) ⋅ forward_jacobian(αβ,p)
    sqrt(g[1,1]*g[2,2] - g[1,2]*g[2,1] )
  end
end

function inv_metric(p)
  function _inv_metric(αβ)
    inv( transpose(forward_jacobian(αβ,p)) ⋅ forward_jacobian(αβ,p) )
  end
end
function inv_metric(αβ,p)
  inv( transpose(forward_jacobian(αβ,p)) ⋅ forward_jacobian(αβ,p) )
end
function measure(αβ,p)
  g = transpose(forward_jacobian(αβ,p)) ⋅ forward_jacobian(αβ,p)
  sqrt(g[1,1]*g[2,2] - g[1,2]*g[2,1] )
end
################################################################################
################################################################################
struct ForwardMap{A}  <: Field
  p::A
end


function Gridap.Arrays.return_cache(f::ForwardMap,cellx::AbstractArray{<:VectorValue{2}})
  p = f.p
  x = first(cellx)
  T = typeof(forward_map(x,p))
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::ForwardMap,cellx::AbstractArray{<:VectorValue{2}} )

  y = cache


  p = f.p

  map!(x -> forward_map(x,p), y, cellx)

  return y
end


function Gridap.Arrays.return_cache(f::ForwardMap,x::VectorValue{2})
  p = f.p

  T = typeof(forward_map(x,p))
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::ForwardMap,x::VectorValue{2})
  y = cache

  p = f.p
  y = forward_map(x,p)

  return y
end

### gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMap},
  cellx::AbstractArray{<:VectorValue{2}})

  p = f.object.p
  x = first(cellx)
  T = typeof(transpose(forward_jacobian(x,p)))

  y = similar(cellx,T)
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:ForwardMap},
  cellx::AbstractArray{<:VectorValue{2}})
  cache, = c
  setsize!(cache,size(cellx))

  y = cache
  p = f.object.p
  map!(x -> transpose(forward_jacobian(x,p)), y, cellx)
  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMap},x::VectorValue{2})
  transpose( forward_jacobian(x,p) )
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:ForwardMap},x::VectorValue{2})
  y = cache
  p = f.object.p
  y = transpose(forward_jacobian(x,p))
  return y
end
