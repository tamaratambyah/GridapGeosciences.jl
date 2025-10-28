const RADIUS = 1.0
function forward_map_2D(p::Int,αβ)
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

forward_map_2D(p::Int) = αβ -> forward_map_2D(p,αβ)

forward_jacobian_2D(p::Int,αβ) = transpose( gradient(forward_map_2D(p))(αβ) )
forward_jacobian_2D(p::Int) = αβ -> forward_jacobian_2D(p,αβ)
