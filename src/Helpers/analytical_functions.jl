################################################################################
#### analytic functions relating to meas --> for debugging purposes
################################################################################
function rho(αβ)
  α,β = αβ
  sqrt(1 + (tan(α))^2 + (tan(β))^2 )
end

function rho3(αβ)
  α,β = αβ
  ( sqrt(1 + (tan(α))^2 + (tan(β))^2 ) )^3
end

function dXda(αβ)
  α,β = αβ
  RADIUS*( - tan(α)*(sec(α))^2 / ( rho3(αβ) ) )
end

function dXdb(αβ)
  α,β = αβ
  RADIUS*( - tan(β)*(sec(β))^2 /(  rho3(αβ)  ) )
end

function dYda(αβ)
  α,β = αβ
  dXda(αβ)*tan(α) + RADIUS/rho(αβ)*(sec(α))^2
end

function dYdb(αβ)
  α,β = αβ
  dXdb(αβ)*tan(α)
end

function dZda(αβ)
  α,β = αβ
  dXda(αβ)*tan(β)
end

function dZdb(αβ)
  α,β = αβ
  dXdb(αβ)*tan(β) + RADIUS/rho(αβ)*(sec(β))^2
end

E(αβ) = dXda(αβ)*dXda(αβ) + dYda(αβ)*dYda(αβ) + dZda(αβ)*dZda(αβ)
F(αβ) = dXda(αβ)*dXdb(αβ) + dYda(αβ)*dYdb(αβ) + dZda(αβ)*dZdb(αβ)
G(αβ) = dXdb(αβ)*dXdb(αβ) + dYdb(αβ)*dYdb(αβ) + dZdb(αβ)*dZdb(αβ)


_sqrtg(αβ) = sqrt( E(αβ)*G(αβ) - F(αβ)*F(αβ) )
