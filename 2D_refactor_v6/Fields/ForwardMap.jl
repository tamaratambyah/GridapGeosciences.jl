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
  drho_da = - tan(α)*(sec(α))^2 / ( rho^3 )
  drho_db = - tan(β)*(sec(β))^2 / ( rho^3 )

  if p == 1
    # X = 1/rho
    # Y = 1/rho * tan(α)
    # Z = 1/rho * tan(β)

    dXda = drho_da
    dXdb = drho_db
    dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
    dYdb = drho_db*tan(α)
    dZda = drho_da*tan(β)
    dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2
  elseif p == 2
    # X = -1/rho * tan(β)
    # Y = 1/rho * tan(α)
    # Z = 1/rho

    dXda = -( drho_da*tan(β) )
    dXdb = -( drho_db*tan(β) + 1/rho*(sec(β))^2 )
    dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
    dYdb = drho_db*tan(α)
    dZda = drho_da
    dZdb = drho_db
  elseif p == 3
    # X = -1/rho * tan(α)
    # Y = 1/rho
    # Z = 1/rho * tan(β)

    dXda = -( drho_da*tan(α) + 1/rho*(sec(α))^2  )
    dXdb = -( drho_db*tan(α) )
    dYda = drho_da
    dYdb = drho_db
    dZda = drho_da*tan(β)
    dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2
  elseif p == 4
    # X = -1/rho
    # Y = 1/rho * tan(β)
    # Z = 1/rho * tan(α)

    dXda = -drho_da
    dXdb = -drho_db
    dYda = drho_da*tan(β)
    dYdb = drho_db*tan(β) + 1/rho*(sec(β))^2
    dZda = drho_da*tan(α) + 1/rho*(sec(α))^2
    dZdb = drho_db*tan(α)
  elseif p == 5
    # X = -1/rho * tan(α)
    # Y = 1/rho * tan(β)
    # Z = -1/rho

    dXda = -(drho_da*tan(α) + 1/rho*(sec(α))^2)
    dXdb = -(drho_db*tan(α) )
    dYda = drho_da*tan(β)
    dYdb = drho_db*tan(β) + 1/rho*(sec(β))^2
    dZda = -drho_da
    dZdb = -drho_db
  elseif p == 6
    # X = -1/rho * tan(β)
    # Y = -1/rho
    # Z = 1/rho * tan(α)

    dXda = -(drho_da*tan(β))
    dXdb = -(drho_db*tan(β) + 1/rho*(sec(β))^2)
    dYda = -drho_da
    dYdb = -drho_db
    dZda = drho_da*tan(α) + 1/rho*(sec(α))^2
    dZdb = drho_db*tan(α)

  end

  ## J = [dXda dXdb
  ##      dYda dYdb
  ##      dZda dZdb  ]
  ## As a TensorValue data = (dXda,dYda,dZda, dXdb, dYdb, dZdb)

  RADIUS*TensorValue{3,2}(dXda,dYda,dZda, dXdb,dYdb,dZdb)

end


## returns 3x2 jacobians
forward_jacobian(p) = αβ -> forward_jacobian(αβ,p)
covarient_basis(p) = αβ -> forward_jacobian(αβ,p)

# return the left psudeo inverse
import Gridap.TensorValues: MultiValue
function pinvJ(J::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  Jt = transpose(J)
  inv(Jt⋅J)⋅Jt
end
forward_pinv_jacobian(p) = αβ -> pinvJ(forward_jacobian(p)(αβ))


# function forward_pinv_jacobian(αβ,p)

#   α,β = αβ
#   rho = sqrt(1 + (tan(α))^2 + (tan(β))^2 )
#   drho_da = - tan(α)*(sec(α))^2 / ( rho^3 )
#   drho_db = - tan(β)*(sec(β))^2 / ( rho^3 )

#   if p == 1
#     dXda = drho_da
#     dXdb = drho_db
#     dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
#     dYdb = drho_db*tan(α)
#     dZda = drho_da*tan(β)
#     dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2
#   elseif p == 2
#     dXda = -( drho_da*tan(β) )
#     dXdb = -( drho_db*tan(β) + 1/rho*(sec(β))^2 )
#     dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
#     dYdb = drho_db*tan(α)
#     dZda = drho_da
#     dZdb = drho_db
#   elseif p == 3
#     dXda = -( drho_da*tan(α) + 1/rho*(sec(α))^2  )
#     dXdb = -( drho_db*tan(α) )
#     dYda = drho_da
#     dYdb = drho_db
#     dZda = drho_da*tan(β)
#     dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2
#   elseif p == 4
#     dXda = -drho_da
#     dXdb = -drho_db
#     dYda = -(-drho_da*tan(β) )
#     dYdb = -(-drho_db*tan(β) + -1/rho*(sec(β))^2 )
#     dZda = -(-drho_da*tan(α) + -1/rho*(sec(α))^2  )
#     dZdb = -(-drho_db*tan(α) )
#   elseif p == 5
#     dXda = -drho_da*tan(α) + -1/rho*(sec(α))^2
#     dXdb = -drho_db*tan(α)
#     dYda = -(-drho_da*tan(β))
#     dYdb = -(-drho_db*tan(β) + -1/rho*(sec(β))^2 )
#     dZda = -drho_da
#     dZdb = -drho_db
#   elseif p == 6
#     dYda = -drho_da
#     dYdb = -drho_db
#     dZda = -(-drho_da*tan(α) + -1/rho*(sec(α))^2  )
#     dZdb = -(-drho_db*tan(α) )
#     dXda = -drho_da*tan(β)
#     dXdb = -drho_db*tan(β) + -1/rho*(sec(β))^2
#   end

#   ## J = [dXda dXdb
#   ##      dYda dYdb
#   ##      dZda dZdb  ]
#   ## As a TensorValue data = (dXda,dYda,dZda, dXdb, dYdb, dZdb)

#   J = RADIUS*TensorValue{3,2}(dXda,dYda,dZda, dXdb,dYdb,dZdb)
#   Jt = RADIUS*TensorValue{2,3}(dXda,dXdb, dYda,dYdb, dZda,dZdb)
#   inv(Jt⋅J)⋅Jt
#   # analytic_inv_metric(αβ) ⋅ Jt
# end

################################################################################
################################################################################
struct ForwardMapPanel1  <: Field
end

"""
The forward map goes 2D -> 3D.
  X = RADIUS/sqrt(1 + (tan(α))^2 + (tan(β))^2  )
  Y = RADIUS/sqrt(1 + (tan(α))^2 + (tan(β))^2  ) * tan(α)
  Z = RADIUS/sqrt(1 + (tan(α))^2 + (tan(β))^2  ) * tan(β)
"""

function Gridap.Arrays.return_cache(f::ForwardMapPanel1,cellx::AbstractArray{<:VectorValue{2}})
  y = similar(cellx,VectorValue{3,Float64})
  return y
end


function Gridap.Arrays.evaluate!(cache,f::ForwardMapPanel1,cellx::AbstractArray{<:VectorValue{2}} )

  y = cache

  map!(x-> Point(RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )),
                 RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )) * tan(x[1]),
                 RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )) * tan(x[2])  ),
      y, cellx  )

  return y
end


function Gridap.Arrays.return_cache(f::ForwardMapPanel1,x::VectorValue{2})
  y = zero(VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::ForwardMapPanel1,x::VectorValue{2})
  y = cache
  y = Point(RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )),
            RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )) * tan(x[1]),
            RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )) * tan(x[2])  )

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
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMapPanel1},
  cellx::AbstractArray{<:VectorValue{2}})
  y = similar(cellx,TensorValue{2,3,Float64})
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:ForwardMapPanel1},
  cellx::AbstractArray{<:VectorValue{2}})
  cache, = c
  setsize!(cache,size(cellx))

  y = cache
  map!(x-> TensorValue{2,3}( dXda(x),dXdb(x),  dYda(x),dYdb(x),  dZda(x),dZdb(x) ),
      y, cellx  )
  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMapPanel1},x::VectorValue{2})
  zero(TensorValue{2,3,Float64})
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:ForwardMapPanel1},x::VectorValue{2})
  y = cache

  y = TensorValue{2,3}( dXda(x),dXdb(x),  dYda(x),dYdb(x),  dZda(x),dZdb(x) )

  return y
end
