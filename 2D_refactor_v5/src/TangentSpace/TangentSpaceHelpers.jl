"""
tangent_f

Compute tangent component of analytical function in ambient space of the sphere
This implementation assumes the normal is (x,y,z)
"""

function tangent_f(X)
  # vec = VectorValue(-X[2],X[1],0)
  vec = VectorValue(X[1]^2,X[1]*X[2],X[3])

  normal_vec = 1/sqrt(X[1]*X[1] + X[2]*X[2] + X[3]*X[3])*VectorValue(X[1],X[2],X[3])
  normal_comp = (vec⋅normal_vec)*normal_vec
  tangent_comp = vec - normal_comp

  @assert norm(normal_vec) ≈ 1.0 # check length of normal vector = 1
  @assert dot(normal_comp, tangent_comp) <= 9.9e-15 # check normal and tangent components are perpendicular

  tangent_comp

end

function tangent_f(vec_func::Function)
  function _tangent_f(X)
    vec = vec_func(X) # VectorValue(-X[2],X[1],0)
    # vec = VectorValue(X[1]^2,X[1]*X[2],X[3])

    normal_vec = 1/sqrt(X[1]*X[1] + X[2]*X[2] + X[3]*X[3])*VectorValue(X[1],X[2],X[3])
    normal_comp = (vec⋅normal_vec)*normal_vec
    tangent_comp = vec - normal_comp

    @assert norm(normal_vec) ≈ 1.0 # check length of normal vector = 1
    @assert dot(normal_comp, tangent_comp) <= 9.9e-15 # check normal and tangent components are perpendicular

    tangent_comp
  end
end


"""
u_parametric

Map analytical vector valued function from X,Y,Z -> α,β
1. Map a point α,β -> X,Y,Z
2. Compute u(X,Y,Z)
3. Compute J cooresponding to α,β -> X,Y,Z
4. Evaluate v(α,β) = pinv(J)⋅u(X,Y,Z)

Note, this is a coordinate transform of an analytic function, so do not need to
apply the piola transform (I think).
Note, the inverse map X,Y,Z -> α,β is ill defined (det(J) = 0). This is why we
need to use the (left) pseudo-inverse of the forward map
"""
function u_parametric(p::Int,uX::Function)
  function _u(αβ)
    cmap = PanelRotationField(r1p_3D[p]) ∘ SigmaField(r) ∘ GnomonicField()

    XYZ = cmap(αβ)

    Jt = ∇(cmap)
    J = Operation(transpose)(Jt) # 3 x 2
    pinvJ = Operation(pinv)(J) # 2 x 3

    Jt_x = Jt(αβ)
    J_x = J(αβ)
    pinvJ_x = pinvJ(αβ)

    # (2x3) (3x1) = (2x1) ∈ α,β
    return pinvJ_x ⋅ uX(XYZ)


  end
end


"""
u_projection

Project an analytical vector valued function from X,Y,Z -> sphere surface
1. Map a point X,Y,Z -> α,β using inverse map
2. Compute u(X,Y,Z)
3. Compute J, J^† cooresponding to α,β -> X,Y,Z
4. Evaluate proj(u) = J ⋅ pinv(J)⋅u(X,Y,Z)

Note, this is a coordinate transform of an analytic function, so do not need to
apply the piola transform (I think)
Note, the inverse map X,Y,Z -> α,β is ill defined for mapping vector fields
det(J) = 0). This is why we need to use the (left) pseudo-inverse of the forward
map. However, it should be okay to map a single point (hopefully)
"""
function u_projection(p::Int,uX::Function)
  function _u(XYZ)
    cmap = PanelRotationField(r1p_3D[p]) ∘ SigmaField(r) ∘ GnomonicField()
    inv_cmap = InvGnomonicField() ∘ InvSigmaField(r)  ∘ PanelRotationField(rp1_3D[p])

    αβ = inv_cmap(XYZ)

    Jt = ∇(cmap)
    J = Operation(transpose)(Jt) # 3 x 2
    pinvJ = Operation(pinv)(J) # 2 x 3

    Jt_x = Jt(αβ)
    J_x = J(αβ)
    pinvJ_x = pinvJ(αβ)

    J_x ⋅ pinvJ_x ⋅ uX(XYZ)

  end
end
