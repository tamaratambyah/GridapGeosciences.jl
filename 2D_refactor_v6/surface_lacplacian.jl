

dfda(f::Function,p::Int) = αβ -> (gradient(f(p))(αβ))[1]
dfdb(f::Function,p::Int) = αβ -> (gradient(f(p))(αβ))[2]

w1(f::Function,p::Int) = αβ -> 1/sqrtg(αβ) * G(αβ)*dfda(f,p)(αβ) - 1/sqrtg(αβ)*F(αβ)*dfdb(f,p)(αβ)
w2(f::Function,p::Int) = αβ -> -1/sqrtg(αβ) * F(αβ)*dfda(f,p)(αβ) + 1/sqrtg(αβ)*E(αβ)*dfdb(f,p)(αβ)

w(f::Function,p::Int) = αβ -> VectorValue(w1(f,p)(αβ),w2(f,p)(αβ))

surflap(f::Function,p::Int) = αβ -> 1/sqrtg(αβ)*(divergence(w(f,p))(αβ))
surflap(f::Function) = p -> surflap(f,p)
