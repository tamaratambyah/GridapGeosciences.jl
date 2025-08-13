### f = sin(ϕ) = Z
function f_sin(p)
  function _f(αβ)
    α,β = αβ
    if p == 1 || p == 3
      return 1/rho(αβ)*tan(β)
    elseif p == 2
      return 1/rho(αβ)
    elseif p == 4 || p == 6
      return 1/rho(αβ)*tan(α)
    elseif p == 5
      return -1/rho(αβ)
    end
  end
end

function rho(αβ)
  α,β = αβ
  sqrt(1 + (tan(α))^2 + (tan(β))^2 )
end

function rho3(αβ)
  α,β = αβ
  ( sqrt(1 + (tan(α))^2 + (tan(β))^2 ) )^3
end



## f = XYZ
function f_XYZ(p)
  function _f(αβ)
    α,β = αβ
    if p == 1 || p == 5 || p == 6
      return RADIUS^3/rho3(αβ)*tan(α)*tan(β)
    else
      return -RADIUS^3/rho3(αβ)*tan(α)*tan(β)
    end
  end
end


################################################################################
#### analytic functions relating to metric, surflap
################################################################################
function dXda(αβ)
  α,β = αβ
  - tan(α)*(sec(α))^2 / ( rho3(αβ)  )
end

function dXdb(αβ)
  α,β = αβ
  - tan(β)*(sec(β))^2 /(  rho3(αβ)  )
end

function dYda(αβ)
  α,β = αβ
  dXda(αβ)*tan(α) + 1/rho(αβ)*(sec(α))^2
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
  dXdb(αβ)*tan(β) + 1/rho(αβ)*(sec(β))^2
end

E(αβ) = dXda(αβ)*dXda(αβ) + dYda(αβ)*dYda(αβ) + dZda(αβ)*dZda(αβ)
F(αβ) = dXda(αβ)*dXdb(αβ) + dYda(αβ)*dYdb(αβ) + dZda(αβ)*dZdb(αβ)
G(αβ) = dXdb(αβ)*dXdb(αβ) + dYdb(αβ)*dYdb(αβ) + dZdb(αβ)*dZdb(αβ)

detg(αβ) = E(αβ)*G(αβ) - F(αβ)*F(αβ)
sqrtg(αβ) = sqrt( E(αβ)*G(αβ) - F(αβ)*F(αβ) )
analytic_metric(αβ) = TensorValue{2,2}(E(αβ),F(αβ),F(αβ),G(αβ))
analytic_inv_metric(αβ) =  TensorValue{2,2}(G(αβ)/detg(αβ),-F(αβ)/detg(αβ),-F(αβ)/detg(αβ),E(αβ)/detg(αβ))
