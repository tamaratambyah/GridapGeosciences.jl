
"""
equi-angular gnomonic map to lat-lon
  input is local cartesian coords on reference panel
  returns lat-lon
"""
function gamma(X::VectorValue{2,Float64})
  x,y = X # local 2D Cartesian coordinates on the reference panel

  @assert x∈[-1.0,1.0] && y∈[-1.0,1.0]

  θ = atan(x)

  bt = (1.0 + x*x + y*y)^(0.5)
  ϕ = asin( y/bt   )

  Point(θ,ϕ)
end
