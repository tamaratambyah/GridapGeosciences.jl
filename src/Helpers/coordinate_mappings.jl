panel_to_cartesian(fX::Function,m::Field) = x -> fX(m(x))
panel_to_cartesian(fX::Function) = m -> panel_to_cartesian(fX,m)

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
