## panel -> XYZ
X_p1(αβ) = RADIUS/rho(αβ)
Y_p1(αβ) = RADIUS/rho(αβ)*tan(αβ[1])
Z_p1(αβ) = RADIUS/rho(αβ)*tan(αβ[2])

#forward_map(p) = αβ -> R1p[p] ⋅ VectorValue(X_p1(αβ),Y_p1(αβ),Z_p1(αβ))

panel_to_cartesian(fX::Function,p::Int) = αβ -> fX(forward_map(p)(αβ))
panel_to_cartesian(fX::Function) = p -> panel_to_cartesian(fX,p)

## panel -> lalon
function latlon_to_cartesian(θϕ::VectorValue{2})
  θ,ϕ = θϕ
  X = cos(θ)*cos(ϕ)
  Y = sin(θ)*cos(ϕ)
  Z = sin(ϕ)
  RADIUS*Point(X,Y,Z)
end

function latlon_to_cartesian_jacobian(θϕ)
  θ,ϕ = θϕ
  dXdtheta = -RADIUS*sin(θ)*cos(ϕ)
  dXdphi = -RADIUS*cos(θ)*sin(ϕ)
  dYdtheta = RADIUS*cos(θ)*cos(ϕ)
  dYdphi = -RADIUS*sin(θ)*sin(ϕ)
  dZdtheta = 0.0
  dZdphi = RADIUS*cos(ϕ)

  ## J = [dXdtheta dXdphi
  ##      dYdtheta dYdphi
  ##      dZdtheta dZdphi  ]
  ## As a TensorValue data = (dXdtheta,dYdtheta,dZdtheta,  dXdphi,dYdphi,dZdphi)
  TensorValue{3,2}(dXdtheta,dYdtheta,dZdtheta,  dXdphi,dYdphi,dZdphi)
end

Jθϕ(XYZ) = latlon_to_cartesian_jacobian(cartesian_to_latlon(XYZ))


function cartesian_to_latlon(XYZ::VectorValue{3})
  X,Y,Z = XYZ
  θ = atan(Y,X)
  ϕ = atan(Z,sqrt(X^2 + Y^2))
  Point(θ,ϕ)
end
cartesian_to_latlon(fθϕ::Function) = XYZ -> fθϕ(cartesian_to_latlon(XYZ))
vec_cartesian_to_latlon(vecθϕ::Function) = XYZ -> Jθϕ(XYZ) ⋅ cartesian_to_latlon(vecθϕ)(XYZ)

panel_to_latlon(f::Function,p::Int) = αβ -> f(cartesian_to_latlon(forward_map(p)(αβ)))
panel_to_latlon(f::Function) = p -> panel_to_latlon(f,p)
