function map_cube_to_sphere(XYZ)
  R = 1
  x,y,z = XYZ
  x_sphere = x*sqrt(1.0-y^2/2-z^2/2+y^2*z^2/(3.0))
  y_sphere = y*sqrt(1.0-z^2/2-x^2/2+x^2*z^2/(3.0))
  z_sphere = z*sqrt(1.0-x^2/2-y^2/2+x^2*y^2/(3.0))
  R*Point(x_sphere,y_sphere,z_sphere)
end

function _map_cube_to_sphere(XYZ)
  R = 1

  x,y,z = XYZ
  factor = R/sqrt(x^2 + y^2 + z^2)
  XYZ*factor
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

function θϕ2xyz(θϕ)
  θ,ϕ = θϕ
  x = cos(θ)*cos(ϕ)
  y = sin(θ)*cos(ϕ)
  z = sin(ϕ)
  VectorValue(x,y,z)
end


function map_unrotated_XYZ(XYZrot,λϕ0)
  λ0,ϕ0 = λϕ0
  Xrot,Yrot,Zrot = XYZrot

  _X = cos(λ0)*sin(ϕ0)*Xrot - sin(λ0)*Yrot + cos(λ0)*cos(ϕ0)*Zrot
  _Y = sin(λ0)*sin(ϕ0)*Xrot + cos(λ0)*Yrot + sin(λ0)*cos(ϕ0)*Zrot
  _Z = - cos(ϕ0)*Xrot + sin(ϕ0)*Zrot
  Point(_X,_Y,_Z)
end

function map_cube_to_latlon(XYZ_cube)
  λϕ0 = Point(0.0,0.0)
  XYZ_sphere = _map_cube_to_sphere(XYZ_cube)
  # _XYZ = map_unrotated_XYZ(XYZ_sphere,λϕ0)
  # xyz2θϕ(_XYZ)
  XYZ_sphere
end
