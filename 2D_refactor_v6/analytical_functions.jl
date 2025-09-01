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

################################################################################
#### analytic functions relating to metric, surflap
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

detg(αβ) = E(αβ)*G(αβ) - F(αβ)*F(αβ)
sqrtg(αβ) = sqrt( E(αβ)*G(αβ) - F(αβ)*F(αβ) )
grad_meas(αβ) = gradient(sqrtg)(αβ)
analytic_metric(αβ) = TensorValue{2,2}(E(αβ),F(αβ),F(αβ),G(αβ))
analytic_inv_metric(αβ) =  TensorValue{2,2}(G(αβ)/detg(αβ),-F(αβ)/detg(αβ),-F(αβ)/detg(αβ),E(αβ)/detg(αβ))

analytic_J1(αβ) = RADIUS*TensorValue{3,2}(dXda(αβ),dYda(αβ),dZda(αβ), dXdb(αβ),dYdb(αβ),dZdb(αβ))

## A = [g12 g22; -g11 -g21] = [F G; -E - F]
## as a TensorValue, (F,-E,G,-F)
analytic_perp_matrix(αβ) = TensorValue{2,2}( F(αβ), -E(αβ), G(αβ), -F(αβ) )

### to compute surflap in components
# dfda(f::Function,p::Int) = αβ -> (gradient(f(p))(αβ))[1]
# dfdb(f::Function,p::Int) = αβ -> (gradient(f(p))(αβ))[2]

# w1(f::Function,p::Int) = αβ -> 1/sqrtg(αβ) * G(αβ)*dfda(f,p)(αβ) - 1/sqrtg(αβ)*F(αβ)*dfdb(f,p)(αβ)
# w2(f::Function,p::Int) = αβ -> -1/sqrtg(αβ) * F(αβ)*dfda(f,p)(αβ) + 1/sqrtg(αβ)*E(αβ)*dfdb(f,p)(αβ)
# w(f::Function,p::Int) = αβ -> VectorValue(w1(f,p)(αβ),w2(f,p)(αβ))
# surflap(f::Function,p::Int) = αβ -> 1/sqrtg(αβ)*(divergence(w(f,p))(αβ))

### to compute surflap more compactly
W(f,p) = αβ ->  sqrtg(αβ)*( analytic_inv_metric(αβ) ⋅ gradient(f(p))(αβ))
surflap(f::Function,p::Int) = αβ -> 1/sqrtg(αβ) * ( divergence(W(f,p))(αβ) )

surflap(f::Function) = p -> surflap(f,p)


### to compute surface divergence analytically
_sdiv(vec::Function,p) = αβ ->  sqrtg(αβ)*( vec(p)(αβ))
surfdiv(vec::Function,p::Int) = αβ -> 1/sqrtg(αβ) * ( divergence(_sdiv(vec,p))(αβ) )
surfdiv(vec::Function) = p -> surfdiv(vec,p)


##### to compute contra sgrad analytically
contr_gradf(f::Function,p::Int) = αβ -> analytic_inv_metric(αβ) ⋅ gradient(f(p))(αβ)
contr_gradf(f::Function) = p -> contr_gradf(f,p)

sgrad(f::Function,p::Int) =  αβ -> forward_jacobian(αβ,p) ⋅ (analytic_inv_metric(αβ) ⋅ gradient(f(p))(αβ))
sgrad(f::Function) = p -> sgrad(f,p)
