const RADIUS = 1.0
function forward_map(p::Int,αβ)
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

forward_map(p::Int) = αβ -> forward_map(p,αβ)

forward_jacobian(p::Int,αβ) = transpose( gradient(forward_map(p))(αβ) )
forward_jacobian(p::Int) = αβ -> forward_jacobian(p,αβ)
covarient_basis(p::Int) = αβ -> forward_jacobian(p,αβ)

# function forward_jacobian(αβ,p)

#   α,β = αβ
#   rho = sqrt(1 + (tan(α))^2 + (tan(β))^2 )
#   drho_da = - tan(α)*(sec(α))^2 / ( rho^3 )
#   drho_db = - tan(β)*(sec(β))^2 / ( rho^3 )

#   if p == 1
#     # X = 1/rho
#     # Y = 1/rho * tan(α)
#     # Z = 1/rho * tan(β)

#     dXda = drho_da
#     dXdb = drho_db
#     dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
#     dYdb = drho_db*tan(α)
#     dZda = drho_da*tan(β)
#     dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2
#   elseif p == 2
#     # X = -1/rho * tan(β)
#     # Y = 1/rho * tan(α)
#     # Z = 1/rho

#     dXda = -( drho_da*tan(β) )
#     dXdb = -( drho_db*tan(β) + 1/rho*(sec(β))^2 )
#     dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
#     dYdb = drho_db*tan(α)
#     dZda = drho_da
#     dZdb = drho_db
#   elseif p == 3
#     # X = -1/rho * tan(α)
#     # Y = 1/rho
#     # Z = 1/rho * tan(β)

#     dXda = -( drho_da*tan(α) + 1/rho*(sec(α))^2  )
#     dXdb = -( drho_db*tan(α) )
#     dYda = drho_da
#     dYdb = drho_db
#     dZda = drho_da*tan(β)
#     dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2
#   elseif p == 4
#     # X = -1/rho
#     # Y = 1/rho * tan(β)
#     # Z = 1/rho * tan(α)

#     dXda = -drho_da
#     dXdb = -drho_db
#     dYda = drho_da*tan(β)
#     dYdb = drho_db*tan(β) + 1/rho*(sec(β))^2
#     dZda = drho_da*tan(α) + 1/rho*(sec(α))^2
#     dZdb = drho_db*tan(α)
#   elseif p == 5
#     # X = -1/rho * tan(α)
#     # Y = 1/rho * tan(β)
#     # Z = -1/rho

#     dXda = -(drho_da*tan(α) + 1/rho*(sec(α))^2)
#     dXdb = -(drho_db*tan(α) )
#     dYda = drho_da*tan(β)
#     dYdb = drho_db*tan(β) + 1/rho*(sec(β))^2
#     dZda = -drho_da
#     dZdb = -drho_db
#   elseif p == 6
#     # X = -1/rho * tan(β)
#     # Y = -1/rho
#     # Z = 1/rho * tan(α)

#     dXda = -(drho_da*tan(β))
#     dXdb = -(drho_db*tan(β) + 1/rho*(sec(β))^2)
#     dYda = -drho_da
#     dYdb = -drho_db
#     dZda = drho_da*tan(α) + 1/rho*(sec(α))^2
#     dZdb = drho_db*tan(α)

#   end

#   ## J = [dXda dXdb
#   ##      dYda dYdb
#   ##      dZda dZdb  ]
#   ## As a TensorValue data = (dXda,dYda,dZda, dXdb, dYdb, dZdb)

#   RADIUS*TensorValue{3,2}(dXda,dYda,dZda, dXdb,dYdb,dZdb)

# end


## returns 3x2 jacobians
# forward_jacobian(p::Int) = αβ -> forward_jacobian(αβ,p)
# covarient_basis(p::Int) = αβ -> forward_jacobian(αβ,p)

# return the left psudeo inverse
import Gridap.TensorValues: MultiValue
function pinvJ(J::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  Jt = transpose(J)
  inv(Jt⋅J)⋅Jt
end
# forward_pinv_jacobian(p) = αβ -> pinvJ(forward_jacobian(p)(αβ))
forward_pinv_jacobian(p) = αβ -> forward_pinv_jacobian(αβ,p)

function forward_pinv_jacobian(αβ,p)

  α,β = αβ

  X,Y,Z = forward_map(p,αβ)
  radius = sqrt(X^2 + Y^2 + Z^2)

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

  J =  radius*TensorValue{3,2}(dXda,dYda,dZda, dXdb,dYdb,dZdb)
  Jt = radius*TensorValue{2,3}(dXda,dXdb, dYda,dYdb, dZda,dZdb)
  inv(Jt⋅J)⋅Jt
  # analytic_inv_metric(αβ) ⋅ Jt
end
