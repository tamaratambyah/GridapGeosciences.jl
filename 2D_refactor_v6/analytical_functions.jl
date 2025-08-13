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
