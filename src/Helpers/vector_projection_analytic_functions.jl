# vector length
vector_length(u::VectorValue{3}) = sqrt(u[1]*u[1] + u[2]*u[2] + u[3]*u[3])
vector_length(u::VectorValue{2}) = sqrt(u[1]*u[1] + u[2]*u[2])

# unit normal
normal_vec(XYZ) = 1/sqrt(XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2] + XYZ[3]*XYZ[3])*VectorValue(XYZ[1],XYZ[2],XYZ[3])

function normal_vector_from_basis(p)
  function _func(αβ)
    J = forward_jacobian(αβ,p)
    a1 = VectorValue(J[1],J[2],J[3])
    a2 = VectorValue(J[4],J[5],J[6])
    n = cross(a1,a2)

    _n = n*(1/sqrtg(αβ) )
    @check vector_length(_n) ≈ 1.0
    _n
  end
end


#### In the following functions, vecX(XYZ) == 3D vector

# tangent component of aribitary 3D vector vecX
tangent_vec(vecX::Function) = XYZ -> vecX(XYZ) - (vecX(XYZ)⋅normal_vec(XYZ))⋅normal_vec(XYZ)


# contravariat components of 3D vector vecX
contra_v(vecX::Function,p::Int) = αβ -> forward_pinv_jacobian(p)(αβ)⋅ vecX(p)(αβ)
contra_v(vecX::Function) = p -> contra_v(vecX,p)


# extract compoents 1 or 2 of contravariat vector, and construct contravariat components of vec perp
contra_v_comp(vecX::Function,p::Int,comp::Int) = αβ -> (forward_pinv_jacobian(p)(αβ)⋅ vecX(p)(αβ))[comp]
contra_v_comp(vecX::Function,comp::Int) = p -> contra_v_comp(vecX,p,comp)

contra_v_perp(vecX::Function,p::Int) = αβ -> sqrtg(αβ)*(
        analytic_inv_metric(αβ) ⋅ VectorValue( -contra_v_comp(vecX,p,2)(αβ), contra_v_comp(vecX,p,1)(αβ) ) )
contra_v_perp(vecX::Function) = p -> contra_v_perp(vecX,p)


# projection of 3D vector vecX
projection_v(vecX::Function,p::Int) = αβ -> forward_jacobian(p)(αβ) ⋅ contra_v(vecX,p)(αβ)
projection_v(vecX::Function) = p -> projection_v(vecX,p)
