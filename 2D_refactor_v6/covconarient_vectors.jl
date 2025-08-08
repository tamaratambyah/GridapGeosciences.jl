using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity
using LinearAlgebra
global RADIUS = 1.0

function forward_map(αβ,p)
  α,β = αβ
  X,Y,Z = 0.0,0.0,0.0

  rho = sqrt( 1 + (tan(α))^2 + (tan(β))^2  )

  if p == 1
    X = 1/rho
    Y = X * tan(α)
    Z = X * tan(β)
  elseif p == 2
    Z = 1/rho
    Y = Z * tan(α)
    X = Z * tan(β)
  elseif p == 3
    Y = 1/rho
    X = Y * tan(α)
    Z = Y * tan(β)
  elseif p == 4
    X = -1/rho
    Y = X * tan(α)
    Z = -X * tan(β)
  elseif p == 5
    Z = -1/rho
    Y = -Z * tan(α)
    X = Z * tan(β)
  elseif p == 6
    Y = -1/rho
    X = Y * tan(α)
    Z = -Y * tan(β)
  end
  XYZ = RADIUS*Point(X,Y,Z)
  return XYZ

end

function forward_jacobian(αβ,p)
  α,β = αβ
  dXda,dYda,dZda = 0.0,0.0,0.0
  dXdb,dYdb,dZdb = 0.0,0.0,0.0

  rho = sqrt( 1 + (tan(α))^2 + (tan(β))^2  )

  dda = (tan(α)*(sec(α))^2 )/ rho^3
  ddb = (tan(β)*(sec(β))^2 )/ rho^3
  factor = 1/rho

  if p == 1
    dXda = -dda
    dXdb = -ddb
    dYda = -dda*tan(α) + factor*(sec(α))^2
    dYdb = -ddb*tan(α)
    dZda = -dda*tan(β)
    dZdb = -ddb*tan(β) + factor*(sec(β))^2
  elseif p == 2
    dZda = -dda
    dZdb = -ddb
    dYda = -dda*tan(α) + factor*(sec(α))^2
    dYdb = -ddb*tan(α)
    dXda = -dda*tan(β)
    dXdb = -ddb*tan(β) + factor*(sec(β))^2
  elseif p == 3
    dYda = -dda
    dYdb = -ddb
    dXda = -dda*tan(α) + factor*(sec(α))^2
    dXdb = -ddb*tan(α)
    dZda = -dda*tan(β)
    dZdb = -ddb*tan(β) + factor*(sec(β))^2
  elseif p == 4
    dXda = dda
    dXdb = ddb
    dYda = dda*tan(α) + -factor*(sec(α))^2
    dYdb = ddb*tan(α)
    dZda = -( dda*tan(β)  )
    dZdb = -( ddb*tan(β) + -factor*(sec(β))^2 )
  elseif p == 5
    dZda = dda
    dZdb = ddb
    dYda = -( dda*tan(α) + -factor*(sec(α))^2 )
    dYdb = -( ddb*tan(α) )
    dXda = dda*tan(β)
    dXdb = ddb*tan(β) + -factor*(sec(β))^2
  elseif p == 6
    dYda = dda
    dYdb = ddb
    dXda = dda*tan(α) + -factor*(sec(α))^2
    dXdb = ddb*tan(α)
    dZda = -( dda*tan(β)  )
    dZdb = -( ddb*tan(β) + -factor*(sec(β))^2 )
  end

  cov_basis_1 = RADIUS*Point(dXda,dYda,dZda)
  cov_basis_2 = RADIUS*Point(dXdb,dYdb,dZdb)
  return cov_basis_1,cov_basis_2

end


function metric(αβ,p)
  cov_basis_1,cov_basis_2 = forward_jacobian(αβ,p)
  E = cov_basis_1 ⋅ cov_basis_1
  F = cov_basis_1 ⋅ cov_basis_2
  G = cov_basis_2 ⋅ cov_basis_2

  g = TensorValue{2,2,Float64}(E,F,F,G)

  return g
end

function meas_metric(αβ,p)
  cov_basis_1,cov_basis_2 = forward_jacobian(αβ,p)
  E = cov_basis_1 ⋅ cov_basis_1
  F = cov_basis_1 ⋅ cov_basis_2
  G = cov_basis_2 ⋅ cov_basis_2

  sqrt_detg = sqrt(E*G - F*F)

  return sqrt_detg
end

function inv_metric(αβ,p)
  cov_basis_1,cov_basis_2 = forward_jacobian(αβ,p)
  E = cov_basis_1 ⋅ cov_basis_1
  F = cov_basis_1 ⋅ cov_basis_2
  G = cov_basis_2 ⋅ cov_basis_2

  detg = (E*G - F*F)

  ginv = 1/detg * TensorValue{2,2}(G,-F,-F,E)

  return ginv
end





function inverse_map(XYZ,p)
  X,Y,Z = XYZ
  α,β = -10.0,-10.0

  if p == 1
    α = atan(Y,X)
    β = atan(Z,X)
  elseif p == 2
    α = atan(Y,Z)
    β = atan(X,Z)
  elseif p == 3
    α = atan(X,Y)
    β = atan(Z,Y)
  elseif p == 4
    α = atan(-Y,-X)
    β = atan(Z,-X)
  elseif p == 5
    α = atan(Y,-Z)
    β = atan(-X,-Z)
  elseif p == 6
    α = atan(-X,-Y)
    β = atan(Z,-Y)
  end

  αβ = Point(α,β)
  return αβ


end

function inverse_jacobian(XYZ,p)
  X,Y,Z = XYZ
  dadX,dadY,dadZ = 0.0,0.0,0.0
  dbdX,dbdY,dbdZ = 0.0,0.0,0.0

  if p == 1
    dadX = -Y/(X^2 + Y^2)
    dadY = X/(X^2 + Y^2)
    dadZ = 0.0
    dbdX = -Z/(X^2 + Z^2)
    dbdY = 0.0
    dbdZ = X/(Z^2 + X^2)
  elseif p == 2
    dadX = 0.0
    dadY = Z/(Z^2 + Y^2)
    dadZ = -Y/(Z^2 + Y^2)
    dbdX = Z/(X^2 + Z^2)
    dbdY = 0.0
    dbdZ = -X/(Z^2 + X^2)
  elseif p == 3
    dadX = Y/(X^2 + Y^2)
    dadY = -X/(X^2 + Y^2)
    dadZ = 0.0
    dbdX = 0.0
    dbdY = -Z/( Z^2 + Y^2)
    dbdZ = Y/(Z^2 + Y^2)
  elseif p == 4
    dadX = -Y/(X^2 + Y^2)
    dadY = X/(X^2 + Y^2)
    dadZ = 0.0
    dbdX = Z/(X^2 + Z^2)
    dbdY = 0.0
    dbdZ = -X/(X^2+Z^2)
  elseif p == 5
    dadX = 0.0
    dadY = -Z/(Z^2 + Y^2)
    dadZ = Y/(Z^2 + Y^2)
    dbdX = Z/(X^2 + Z^2)
    dbdY = 0.0
    dbdZ = -X/(Z^2 + X^2)
  elseif p == 6
    dadX = Y/(X^2 + Y^2)
    dadY = -X/(X^2 + Y^2)
    dadZ = 0.0
    dbdX = 0.0
    dbdY = Z/( Z^2 + Y^2)
    dbdZ = -Y/(Z^2 + Y^2)
  end

  contra_basis_1 = Point(dadX,dadY,dadZ)
  contra_basis_2 = Point(dbdX,dbdY,dbdZ)
  return contra_basis_1,contra_basis_2
end

function J(p)
  function _J(αβ)
    XYZ = forward_map(αβ,p)
    contra_basis_1,contra_basis_2 = inverse_jacobian(XYZ,p)
    TensorValue{3,2}(contra_basis_1[1],contra_basis_1[2],contra_basis_1[3],
                    contra_basis_2[1],contra_basis_2[2],contra_basis_2[3]  )
  end
end

function Jcov(p)
  function _J(αβ)
    cov_basis_1,cov_basis_2 = forward_jacobian(αβ,p)
    TensorValue{3,2}(cov_basis_1[1],cov_basis_1[2],cov_basis_1[3],
                    cov_basis_2[1],cov_basis_2[2],cov_basis_2[3]  )
  end
end

function a1(p)
  function _a1(αβ)
    XYZ = forward_map(αβ,p)
    contra_basis_1,contra_basis_2 = inverse_jacobian(XYZ,p)
    return contra_basis_1
  end
end

function convarient_components(vec_phys,αβ,p)
  cov_basis_1,cov_basis_2 = forward_jacobian(αβ,p)
  u_cov = vec_phys ⋅ cov_basis_1
  v_cov = vec_phys ⋅ cov_basis_2
  return u_cov,v_cov
end

function contrvarient_components(vec_phys,XYZ,p)
  contra_basis_1,contra_basis_2 = inverse_jacobian(XYZ,p)
  u_contra = vec_phys ⋅ contra_basis_1
  v_contra = vec_phys ⋅ contra_basis_2
  return u_contra,v_contra
end


function cov_to_phys(vec_phys,αβ,p)
  XYZ = forward_map(αβ,p)
  u_cov,v_cov = convarient_components(vec_phys,αβ,p)

  contra_basis_1,contra_basis_2 = inverse_jacobian(XYZ,p)
  u = u_cov*contra_basis_1[1] + v_cov*contra_basis_2[1]
  v = u_cov*contra_basis_1[2] + v_cov*contra_basis_2[2]
  w = u_cov*contra_basis_1[3] + v_cov*contra_basis_2[3]

  return Point(u,v,w)

end

function contra_to_phys(vec_phys,αβ,p)
  XYZ = forward_map(αβ,p)
  u_contra,v_contra = contrvarient_components(vec_phys,XYZ,p)
  cov_basis_1,cov_basis_2 = forward_jacobian(αβ,p)

  u = u_contra*cov_basis_1[1] + v_contra*cov_basis_2[1]
  v = u_contra*cov_basis_1[2] + v_contra*cov_basis_2[2]
  w = u_contra*cov_basis_1[3] + v_contra*cov_basis_2[3]

  return Point(u,v,w)
end


######## mapping of tangent vector. At node 4
vec = Point(-2.0,1.0,1.0)

# panel 1:
p = 1
αβ = Point(π/4,π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 3:
p = 3
αβ = Point(π/4,π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 2:
p = 2
αβ = (π/4,π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

######## mapping of tangent vector. At node 3
vec = Point(1.0,2.0,1.0)

# panel 1:
p = 1
αβ = Point(-π/4,π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 6:
p = 6
αβ = Point(-π/4,π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 2:
p = 2
αβ = (-π/4,π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)


######## mapping of tangent vector. At node 5
vec = Point(1.0,1.0,2.0)

# panel 4:
p = 4
αβ = Point(π/4,π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 6:
p = 6
αβ = Point(π/4,π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

p = 2
αβ = (-π/4,-π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)


######## mapping of tangent vector. At node 1
vec = Point(2.0,1.0,1.0)

# panel 1:
p = 1
αβ = Point(-π/4,-π/4)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 5:
p = 5
αβ = Point(-π/4,-π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)



######## mapping of tangent vector. At node 2
vec = Point(1.0,1.0,2.0)

# panel 1:
p = 1
αβ = Point(π/4,-π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 5:
p = 5
αβ = Point(π/4,-π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)


######## mapping of tangent vector. At node 8
vec = Point(-1.0,-1.0,2.0)

# panel 1:
p = 6
αβ = Point(π/4,-π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 4:
p = 4
αβ = Point(π/4,-π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 5:
p = 5
αβ = Point(-π/4,π/4)
XYZ = forward_map(αβ,p)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)
