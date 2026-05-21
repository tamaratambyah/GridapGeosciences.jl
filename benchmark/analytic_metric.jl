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

inv_sqrtg(αβ) = 1/_sqrtg(αβ)

g(αβ) = TensorValue(E(αβ), F(αβ), F(αβ), G(αβ)  )

inv_g(αβ) = 1/_detg(αβ)* TensorValue(G(αβ), -F(αβ), -F(αβ), E(αβ)  )

function allg(αβ)

  t1 = tan(αβ[1])
  t2 = tan(αβ[2])

  c1 = cos(αβ[1])
  c2 = cos(αβ[2])

  r = 1.0 + t1*t1 + t2*t2
  r4 = r*r

  c = 1.0/( r4*c1*c1 *c2*c2 )

  f = c * ( -1.0*t1*t2 )
  e = c * ( 1.0 + (t1*t1) )
  g = c * ( 1 + (t2*t2) )

  sqrtg = (e*g - f*f)^(1/2)

  1/sqrtg*TensorValue(g, -f, -f, e  )
end


A = [1 0 0
    0 1 0
    0 0 1]
a(x) = 1.0 + 0.0*x[1] + 0.0*x[1]
b(x) = 0.0 + 0.0*x[1] + 0.0*x[1]
id_sqrtg(x) = sqrt( a(x)*a(x) - b(x)*b(x) )
id_detg(x) = a(x)*a(x) - b(x)*b(x)
identity_g(x) = 1/id_detg(x)* TensorValue(a(x),b(x),b(x),  b(x),a(x),b(x),  b(x),b(x),a(x)   )
identity_g(Point(1,1,1))
