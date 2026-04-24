# unit normal
normal_vec(XYZ) = 1/sqrt(XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2] + XYZ[3]*XYZ[3])*VectorValue(XYZ[1],XYZ[2],XYZ[3])



#### In the following functions, vecX(XYZ) == 3D vector

# tangent component of aribitary 3D vector vecX
tangent_vec(vecX::Function) = XYZ -> vecX(XYZ) - (vecX(XYZ)⋅normal_vec(XYZ))⋅normal_vec(XYZ)

# Piola contravariant components of a 3D vector vecX
# The Piola map is ̃u = J ( 1/√g u)
# so u = √g J^† ̃u
piola(vecX::Function,p::Int) = αβ -> sqrtg(p,αβ)*( forward_pinv_jacobian(p)(αβ)⋅ vecX(p)(αβ))
piola(vecX::Function) = p -> piola(vecX,p)

piola(vecX::Function,m::ForwardMap2Dor3D) = αβ -> sqrtg(m,αβ)*( forward_pinv_jacobian(m)(αβ)⋅vecX(m)(αβ))

# Contravariant components of 3D vector vecX
# The contravariatn mapping is  ̃u = J u
# so u = J^† ̃u
contra_v(vecX::Function,p::Int) = αβ -> forward_pinv_jacobian(p)(αβ)⋅vecX(p)(αβ)
contra_v(vecX::Function) = p -> contra_v(vecX,p)

contra_v(vecX::Function,m::ForwardMap2Dor3D) = αβ -> forward_pinv_jacobian(m)(αβ)⋅vecX(m)(αβ)


# projection of 3D vector vecX
projection_v(vecX::Function,p::Int) = αβ -> forward_jacobian(p)(αβ) ⋅ contra_v(vecX,p)(αβ)
projection_v(vecX::Function) = p -> projection_v(vecX,p)
