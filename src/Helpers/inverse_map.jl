function inverse_map(XYZ,p)
  X,Y,Z = XYZ

  if p == 1
    α = atan(Y,X)
    β = atan(Z,X)
  elseif p == 2
    α = atan(Y,Z)
    β = -atan(X,Z)
  elseif p == 3
    α = -atan(X,Y)
    β = atan(Z,Y)
  elseif p == 4
    α = atan(Z,-X)
    β = -atan(-Y,-X)
  elseif p == 5
    α = -atan(X,-Z)
    β = atan(Y,-Z)
  elseif p == 6
    α = atan(Z,-Y)
    β = -atan(X,-Y)
  end
  Point(α,β)

end



function inverse_jacobian(XYZ,p)
  X,Y,Z = XYZ

  if p == 1
    dadX = - Y/(X^2 + Y^2)
    dadY = X/(X^2 + Y^2)
    dadZ = 0.0
    dbdX = -Z/(X^2 + Z^2)
    dbdY = 0.0
    dbdZ = X/(X^2 + Z^2)
  elseif p == 2
    dadX = 0.0
    dadY = Z/(Y^2 + Z^2)
    dadZ = -Y/(Y^2 + Z^2)
    dbdX = -Z/(X^2 + Z^2)
    dbdY = 0.0
    dbdZ = X/(X^2 + Z^2)
  elseif p == 3
    dadX = -Y/(X^2 + Y^2)
    dadY = X/(X^2 + Y^2)
    dadZ = 0.0
    dbdX = 0.0
    dbdY = -Z/(Y^2 + Z^2)
    dbdZ = Y/(Y^2 + Z^2)
  elseif p == 4
    dadX = Z/(Z^2 + X^2)
    dadY = 0.0
    dadZ = -X/(Z^2 + X^2)
    dbdX = Y/(X^2 + Y^2)
    dbdY = -X/(X^2 + Y^2)
    dbdZ = 0.0
  elseif p == 5
    dadX = Z/(X^2 + Z^2)
    dadY = 0.0
    dadZ = -X/(X^2 + Z^2)
    dbdX = 0.0
    dbdY = -Z/(Y^2 + Z^2)
    dbdZ = Y/(Y^2 + Z^2)
  elseif p == 6
    dadX = 0.0
    dadY = Z/(Y^2 + Z^2)
    dadZ = -Y/(Y^2 + Z^2)
    dbdX = Y/(X^2 + Y^2)
    dbdY = -X/(X^2 + Y^2)
    dbdZ = 0.0
  end

  ## J = [dadX dadY dadX
  ##      dbdX dbdY dbdZ ]
  ## As a TensorValue data = (dadX,dadY,dadZ, dbdX, dbdY, dbdZ)

  TensorValue{2,3}(dadX,dbdX, dadY,dbdY, dadZ,dbdZ)

end

## returns 2x3 jacobian
function inverse_jacobian(p)
  function _inverse_jacobian(αβ)
    XYZ = forward_map(αβ,p)
    inverse_jacobian(XYZ,p)
  end
end

## returns 3x2 jacobian
function contravariant_basis(p)
  function _contravariant_basis(αβ)
    XYZ = forward_map(αβ,p)
    transpose(inverse_jacobian(XYZ,p))
  end
end
