
## panel -> XYZ
X_p1(αβ) = RADIUS/rho(αβ)
Y_p1(αβ) = RADIUS/rho(αβ)*tan(αβ[1])
Z_p1(αβ) = RADIUS/rho(αβ)*tan(αβ[2])

forward_map(p) = αβ -> R1p[p] ⋅ VectorValue(X_p1(αβ),Y_p1(αβ),Z_p1(αβ))

panel_to_cartesian(f::Function,p::Int) = αβ -> f(forward_map(p)(αβ))
panel_to_cartesian(f::Function) = p -> panel_to_cartesian(f,p)

## panel -> lalon
function latlon_to_cartesian(θϕ::VectorValue{2})
  θ,ϕ = θϕ
  X = cos(θ)*cos(ϕ)
  Y = sin(θ)*cos(ϕ)
  Z = sin(ϕ)
  RADIUS*Point(X,Y,Z)
end

function cartesian_to_latlon(XYZ::VectorValue{3})
  X,Y,Z = XYZ
  θ = atan(Y,X)
  ϕ = atan(Z,sqrt(X^2 + Y^2))
  Point(θ,ϕ)
end

panel_to_latlon(f::Function,p::Int) = αβ -> f(cartesian_to_latlon(forward_map(p)(αβ)))
panel_to_latlon(f::Function) = p -> panel_to_latlon(f,p)
