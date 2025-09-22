### f = sin(ϕ) = Z
function f_sin(p)
  function _f(αβ)
    α,β = αβ
    if p == 1 || p == 3
      return RADIUS/rho(αβ)*tan(β)
    elseif p == 2
      return RADIUS/rho(αβ)
    elseif p == 4 || p == 6
      return RADIUS/rho(αβ)*tan(α)
    elseif p == 5
      return -RADIUS/rho(αβ)
    end
  end
end

### f = sin(ϕ)
fθϕ(θϕ) = sin(θϕ[2])

### f = streamfunction
fWilliamson(ζ) = θϕ -> - sin(θϕ[2])*cos(ζ) + cos(θϕ[1])*cos(θϕ[2])*sin(ζ)


### f = XYZ
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

fX(XYZ::VectorValue{3}) = XYZ[1]*XYZ[2]*XYZ[3]
