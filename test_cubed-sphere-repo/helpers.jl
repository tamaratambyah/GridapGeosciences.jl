# projection of 3D vector vecX
projection_v(vecX::Function,p::Field) = αβ -> forward_jacobian(p)(αβ) ⋅ contra_v(vecX,p)(αβ)
projection_v(vecX::Function) = p -> projection_v(vecX,p)
