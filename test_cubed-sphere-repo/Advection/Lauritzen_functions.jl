"""
initial conditions from Lauritzen2012 paper
doi:10.5194/gmd-5-887-2012
"""
RADIUS = 1.0
θϕc1 = Point(5*π/6,0)
θϕc2 = Point(7*π/6,0)
radius_func(θϕ,θϕc) = RADIUS*acos( sin(θϕc[2])*sin(θϕ[2]) + cos(θϕc[2])*cos(θϕ[2])*cos(θϕ[1]-θϕc[1])   )

### Scalar fields: the exact solution is ϕ(t=T) = ϕ(t=0) where T = 5
function gaussian_hill(xyz)
  h, b = 0.95, 5.0
  x,y,z = xyz
  x1,y1,z1 = θϕ2xyz(θϕc1)
  x2,y2,z2 = θϕ2xyz(θϕc2)

  h*exp( -b*( (x-x1)^2 + (y-y1)^2 + (z-z1)^2 ) ) + h*exp( -b*( (x-x2)^2 + (y-y2)^2 + (z-z2)^2 ) )
end

function cosine_bell(xyz)
  b,c = 0.1, 0.9

  X,Y,Z = xyz
  radius = sqrt(X^2 + Y^2 + Z^2)
  r = radius/2

  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,_r = θϕr

  r1 = radius_func(Point(θ,ϕ),θϕc1)
  r2 = radius_func(Point(θ,ϕ),θϕc2)

  if r1 < r && r2 < r
    return b
  elseif r1 < r
    return b + c*( 0.5*(1 + cos(π*r1/r) ) )
  elseif r2 < r
    return b + c*( 0.5*(1 + cos(π*r2/r) ) )
  else
    return b
  end

end

function slotted_cylinders(xyz)
  b,c = 0.1, 1.0

  X,Y,Z = xyz
  radius = sqrt(X^2 + Y^2 + Z^2)
  r = radius/2

  θ1,ϕ1 = θϕc1
  θ2,ϕ2 = θϕc2

  θ,ϕ,_r   = xyz2θϕr(xyz)

  r1 = radius_func(Point(θ,ϕ),θϕc1)
  r2 = radius_func(Point(θ,ϕ),θϕc2)

  if (r1 <= r) && (abs(θ-θ1) >= r/(6*radius) )
    return c
  elseif (r2 <= r) && (abs(θ-θ2) >= r/(6*radius) )
    return c
  elseif (r1 <= r) && (abs(θ-θ1) < r/(6*radius) ) && ( (ϕ-ϕ1) < -5/12*r/radius)
    return c
  elseif (r2 <= r) && (abs(θ-θ2) < r/(6*radius) ) && ( (ϕ-ϕ2) > 5/12*r/radius)
    return c
  else
    return b
  end

end

function cosine_bell_squared(xyz)
  b,c = 0.1, 0.9

  X,Y,Z = xyz
  radius = sqrt(X^2 + Y^2 + Z^2)
  r = radius/2

  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,_r = θϕr

  r1 = radius_func(Point(θ,ϕ),θϕc1)
  r2 = radius_func(Point(θ,ϕ),θϕc2)

  if r1 < r && r2 < r
    return (b)^2
  elseif r1 < r
    return (b + c*( 0.5*(1 + cos(π*r1/r) ) ))^2
  elseif r2 < r
    return (b + c*( 0.5*(1 + cos(π*r2/r) ) ))^2
  else
    return (b)^2
  end

end

function correlated_cosine_bell(xyz)
  a, b = -0.8, 0.9
  a*(cosine_bell(xyz))^2 + b
  # a*cosine_bell_squared(xyz) + b
end

### Velocity fields
function nondivergent_velocity(t::Float64)
  function _velocity(xyz)
    R, T = 1.0, 5.0
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,_r = θϕr
    u     = 10*R/T*( sin(θ-2*π*t/T) )^2*sin(2*ϕ)*cos(π*t/T) + 2*π*R/T*cos(ϕ)
    v     = 10*R/T*sin(2*(θ-2*π*t/T))*cos(ϕ)*cos(π*t/T)
    _spherical_to_cartesian_matrix(θϕr)⋅VectorValue(u,v,0)
  end
end

function divergent_velocity(t::Float64)
  function _velocity(xyz)
    R, T = 1.0, 5.0
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,_r = θϕr
    u     = -5*R/T*( sin( (θ-2*π*t/T))/2 )^2*sin(2*ϕ)*(cos(ϕ))^2*cos(π*t/T) + 2*π*R/T*cos(ϕ)
    v     = (5/2)*R/T*sin(θ-2*π*t/T)*(cos(ϕ))^3*cos(π*t/T)
    _spherical_to_cartesian_matrix(θϕr)⋅VectorValue(u,v,0)
  end
end

function _spherical_to_cartesian_matrix(θϕr)
  θ,ϕ,r = θϕr
  TensorValue(-sin(θ)       , cos(θ)       ,      0,
              -sin(ϕ)*cos(θ),-sin(ϕ)*sin(θ), cos(ϕ),
               cos(ϕ)*cos(θ), cos(ϕ)*sin(θ), sin(ϕ))
end
