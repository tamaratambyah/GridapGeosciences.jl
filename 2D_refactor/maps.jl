function map_cube_to_sphere(XYZ)
  R = 1
  x,y,z = XYZ
  x_sphere = x*sqrt(1.0-y^2/2-z^2/2+y^2*z^2/(3.0))
  y_sphere = y*sqrt(1.0-z^2/2-x^2/2+x^2*z^2/(3.0))
  z_sphere = z*sqrt(1.0-x^2/2-y^2/2+x^2*y^2/(3.0))
  R*Point(x_sphere,y_sphere,z_sphere)
end


function map_cube_to_cube(XYZ)
  x,y,z = XYZ
  Point(2*x,4*y,-3*z)
end

function xyz2θϕ(x)
  θ = atan(x[2], x[1])
  ϕ = asin(x[3])
  VectorValue(θ,ϕ)
end


function map_cube_to_latlon(XYZ)
  xyz2θϕ(XYZ)
end
