# vector length
vector_length(u::VectorValue{3}) = sqrt(u[1]*u[1] + u[2]*u[2] + u[3]*u[3])
vector_length(u::VectorValue{2}) = sqrt(u[1]*u[1] + u[2]*u[2])

# unit normal
normal_vec(XYZ) = 1/sqrt(XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2] + XYZ[3]*XYZ[3])*VectorValue(XYZ[1],XYZ[2],XYZ[3])

function normal_vector_from_basis(forward_map)
  function _func(αβ)
    Jac = J(forward_map,αβ)
    a1 = VectorValue(Jac[1],Jac[2],Jac[3])
    a2 = VectorValue(Jac[4],Jac[5],Jac[6])
    n = cross(a1,a2)
    _n = n*(1/sqrtg(forward_map,αβ) )
    @check vector_length(_n) ≈ 1.0
    _n
  end
end


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



# extract compoents 1 or 2 of contravariat vector, and construct contravariat components of vec perp
contra_v_comp(vecX::Function,p::Int,comp::Int) = αβ -> (forward_pinv_jacobian(p)(αβ)⋅ vecX(p)(αβ))[comp]
contra_v_comp(vecX::Function,comp::Int) = p -> contra_v_comp(vecX,p,comp)

contra_v_perp(vecX::Function,p::Int) = αβ -> sqrtg(p,αβ)*(
        inv_metric(p,αβ) ⋅ VectorValue( -contra_v_comp(vecX,p,2)(αβ), contra_v_comp(vecX,p,1)(αβ) ) )
contra_v_perp(vecX::Function) = p -> contra_v_perp(vecX,p)


# projection of 3D vector vecX
projection_v(vecX::Function,p::Int) = αβ -> forward_jacobian(p)(αβ) ⋅ contra_v(vecX,p)(αβ)
projection_v(vecX::Function) = p -> projection_v(vecX,p)


#### contravariant components for 3D thick shere
contra_v_comp3D(vecX::Function,p::Int,comp::Int) = αβ -> (forward_pinv_jacobian(p)(αβ)⋅ vecX(p)(αβ))[comp]
contra_v_comp3D(vecX::Function,comp::Int) = p -> contra_v_comp3D(vecX,p,comp)

#### contravariant components of perp for 3D thick shere
contra_v_perp3D(vecX::Function,p::Int) = αβ -> sqrtg(p,αβ)*(
        inv_metric(p,αβ) ⋅ VectorValue(0.0, -contra_v_comp3D(vecX,p,3)(αβ), contra_v_comp3D(vecX,p,2)(αβ) ) )
contra_v_perp3D(vecX::Function) = p -> contra_v_perp3D(vecX,p)
