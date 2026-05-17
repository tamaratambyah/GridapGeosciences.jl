#### Analytic metric
function rho4(αβ)
  α,β = αβ
  (1 + (tan(α))^2 + (tan(β))^2 )^2
end

coeff(αβ) = radius^2/( rho4(αβ) *(cos(αβ[1]))^2 *(cos(αβ[2]))^2 )
E(αβ) = coeff(αβ) * ( 1 + (tan(αβ[1]))^2 )
F(αβ) = coeff(αβ) * ( -1.0*tan(αβ[1])*tan(αβ[2]) )
G(αβ) = coeff(αβ) * ( 1 + (tan(αβ[2]))^2 )

_sqrtg(αβ) = sqrt( E(αβ)*G(αβ) - F(αβ)*F(αβ) )
_detg(αβ) =  E(αβ)*G(αβ) - F(αβ)*F(αβ)

inv_g(αβ) = 1/_detg(αβ)* TensorValue(G(αβ), -F(αβ), -F(αβ), E(αβ)  )
