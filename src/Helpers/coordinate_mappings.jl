panel_to_cartesian(fX::Function,p::Int) = αβ -> fX(forward_map(p)(αβ))
panel_to_cartesian(fX::Function) = p -> panel_to_cartesian(fX,p)

function cartesian_to_latlon(XYZ::VectorValue{3})
  X,Y,Z = XYZ
  θ = atan(Y,X)
  ϕ = atan(Z,sqrt(X^2 + Y^2))
  Point(θ,ϕ)
end
cartesian_to_latlon(fθϕ::Function) = XYZ -> fθϕ(cartesian_to_latlon(XYZ))

panel_to_latlon(f::Function,p::Int) = αβ -> f(cartesian_to_latlon(forward_map(p)(αβ)))
panel_to_latlon(f::Function) = p -> panel_to_latlon(f,p)

################################################################################
### Alberto's coordinate_mappings
################################################################################
function xyz2θϕ(x)
  r = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
  θ = atan(x[2], x[1])
  ϕ = asin(x[3]/r)
  VectorValue(θ,ϕ)
end

function xyz2θϕr(x)
  r = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
  θ = atan(x[2], x[1])
  ϕ = asin(x[3]/r)
  VectorValue(θ,ϕ,r)
end

function θϕ2xyz(θϕ)
  θ,ϕ = θϕ
  x = cos(θ)*cos(ϕ)
  y = sin(θ)*cos(ϕ)
  z = sin(ϕ)
  VectorValue(x,y,z)
end

function spherical_to_cartesian_matrix(θϕr)
  θ,ϕ,r = θϕr
  TensorValue(-sin(θ)       , cos(θ)       ,      0,
              -sin(ϕ)*cos(θ),-sin(ϕ)*sin(θ), cos(ϕ),
               cos(ϕ)*cos(θ), cos(ϕ)*sin(θ), sin(ϕ))
end
