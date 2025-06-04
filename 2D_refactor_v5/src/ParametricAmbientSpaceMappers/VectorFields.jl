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



"""
u_vector_ambient2parametric

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
function u_vector_ambient2parametric(p::Int,uX::Function)
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
u_vector_latlon2ambient

Map analytical vector valued function from θ,ϕ -> X,Y,Z
1. Map a point X,Y,Z -> θ,ϕ
2. Compute u(θ,ϕ)
3. Compute J cooresponding to θ,ϕ -> X,Y,Z (this map is well posed)
4. Evaluate v(X,Y,Z) = J⋅u(θ,ϕ)

Note, this is a coordinate transform of an analytic function, so do not need to
apply the piola transform (I think).
Note, the map X,Y,Z -> θ,ϕ is well defined. This is why we use the J of the
forward map
"""
function u_vector_latlon2ambient(uθϕ::Function)
  function _u(XYZ)
    θϕ = InvSigmaField(r)(XYZ)

    Jt = gradient(SigmaField(r)) # 2 x 3
    J = Operation(transpose)(Jt) # 3 x 2

    J(θϕ) ⋅ uθϕ(θϕ)
  end
end




function parametric_cf_2_ambient_vector(manifold_model,cf_parametric::CellField)
  ambient_model = get_ambient_model(manifold_model)
  panel_ids = get_panel_ids(manifold_model)

  Ω_parametric = Triangulation(manifold_model)
  Ω_ambient = Triangulation(ambient_model)

  mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)
  inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)

  _cf = change_domain(cf_parametric,ReferenceDomain(),PhysicalDomain())

  ## map parametric FEFunction to ambient space
  Jt = lazy_map(Broadcasting(gradient),mapping)
  J = lazy_map(Operation(transpose),Jt)

  _cf_mapped = lazy_map(Broadcasting(⋅),J,get_data(_cf))
  cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,inv_mapping)

  cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() ) # ambient cell field

  return cf_ambient
end
