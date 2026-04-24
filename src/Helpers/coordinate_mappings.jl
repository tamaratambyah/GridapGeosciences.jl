panel_to_cartesian(fX::Function,m::Field) = x -> fX(m(x))
panel_to_cartesian(fX::Function) = m -> panel_to_cartesian(fX,m)

function cartesian_to_latlon(XYZ::VectorValue{3})
  X,Y,Z = XYZ
  θ = atan(Y,X)
  ϕ = atan(Z,sqrt(X^2 + Y^2))
  Point(θ,ϕ)
end
cartesian_to_latlon(fθϕ::Function) = XYZ -> fθϕ(cartesian_to_latlon(XYZ))

panel_to_latlon(f::Function,p::Int) = αβ -> f(cartesian_to_latlon(forward_maps[p](αβ)))
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

### This is Alberto's original spherical_to_cartesian_matrix
### I'm pretty sure it is incorrect.
###   X = r cosθ cosϕ; Y = r sinθ sinϕ ; Z = r sinϕ
###   TensorValue = (  ?,?,dZdθ   dXdϕ,dYdϕ,dZdϕ  , dXdr,dYdr,dZdr)
### Replacing with mine below
# function spherical_to_cartesian_matrix(θϕr)
#   θ,ϕ,r = θϕr
#   TensorValue(-sin(θ)       , cos(θ)       ,      0,
#               -sin(ϕ)*cos(θ),-sin(ϕ)*sin(θ), cos(ϕ),
#                cos(ϕ)*cos(θ), cos(ϕ)*sin(θ), sin(ϕ))
# end

##

"""
spherical_to_cartesian_matrix
  X = r cosθ cosϕ
  Y = r sinθ sinϕ
  Z = r sinϕ
  J = [dXdθ dXdϕ dXdr
       dYdθ dYdϕ dYdr
       dZdθ dZdϕ dZdr ]
As a TensorValue:
  TensorValue = (dXdθ,dYdθ,dZdθ,  dXdϕ,dYdϕ,dZdϕ,  dXdr,dYdr,dZdr)
"""
function spherical_to_cartesian_matrix(θϕr)
  θ,ϕ,r = θϕr
  TensorValue(-r*sin(θ)*cos(ϕ), r*cos(θ)*cos(ϕ), 0,
              -r*cos(θ)*sin(ϕ), -r*sin(θ)*sin(ϕ), r*cos(ϕ),
                 cos(θ)*cos(ϕ),    sin(θ)*cos(ϕ), sin(ϕ) )
end
