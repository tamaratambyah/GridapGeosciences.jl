#### In the following functions, vecX(XYZ) == 3D vector

# tangent component of aribitary 3D vector vecX
tangent_vec(vecX::Function) = XYZ -> vecX(XYZ) - (vecX(XYZ)⋅normal_vec(XYZ))⋅normal_vec(XYZ)

# Piola contravariant components of a 3D vector vecX
# The Piola map is ̃u = J ( 1/√g u)
# so u = √g J^† ̃u
piola(vecX::Function,m::Field) = αβ -> sqrtg(m,αβ)*( forward_pinv_jacobian(m)(αβ)⋅vecX(m)(αβ))
piola(vecX::Function) = p -> piola(vecX,p)

# Contravariant components of 3D vector vecX
# The contravariatn mapping is  ̃u = J u
# so u = J^† ̃u
contra_v(vecX::Function,m::Field) = αβ -> forward_pinv_jacobian(m)(αβ)⋅vecX(m)(αβ)
contra_v(vecX::Function) = p -> contra_v(vecX,p)
