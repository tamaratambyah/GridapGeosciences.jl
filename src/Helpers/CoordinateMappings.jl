panel_to_cartesian(fX::Function,m::Field) = x -> fX(m(x))
panel_to_cartesian(fX::Function) = m -> panel_to_cartesian(fX,m)

################################################################################
### Coordinate_mappings
################################################################################
function xyz2θϕr(x)
  r = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
  θ = atan(x[2], x[1])
  ϕ = asin(x[3]/r)
  VectorValue(θ,ϕ,r)
end
