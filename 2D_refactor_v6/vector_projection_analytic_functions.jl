# unit normal
normal_vec(XYZ) = 1/sqrt(XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2] + XYZ[3]*XYZ[3])*VectorValue(XYZ[1],XYZ[2],XYZ[3])

# vecX(XYZ) == 3D vecotr

# tangent component of aribitary 3D vector vecX
tangent_vec(vecX::Function) = XYZ -> vecX(XYZ) - (vecX(XYZ)⋅normal_vec(XYZ))⋅normal_vec(XYZ)


# contravariat components of 3D vector vecX
contra_v(vecX::Function,p::Int) = αβ -> forward_pinv_jacobian(p)(αβ)⋅ vecX(p)(αβ)
contra_v(vecX::Function) = p -> contra_v(vecX,p)

# projection of 3D vector vecX
projection_v(vecX::Function,p::Int) = αβ -> forward_jacobian(p)(αβ) ⋅ contra_v(vecX,p)(αβ)
projection_v(vecX::Function) = p -> projection_v(vecX,p)
