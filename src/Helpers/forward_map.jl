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


# return the left psudeo inverse
function pinvJ(J::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  Jt = transpose(J)
  inv(Jt⋅J)⋅Jt
end

forward_pinv_jacobian(p) = αβ -> pinvJ(forward_jacobian(p)(αβ))
