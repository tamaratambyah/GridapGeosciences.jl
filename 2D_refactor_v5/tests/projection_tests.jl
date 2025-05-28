using Gridap
include("../src/initialise.jl")
include("../src/Visualization/VisualizationData.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

ambient_model = get_ambient_model(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)

################################################################################
##### Interoplate cell-wise array of generic fields
################################################################################
function Gridap.FESpaces._cell_vals(fs::SingleFieldFESpace,object::AbstractArray{<:GenericField})
  println("my cell vals")
  s = get_fe_dof_basis(fs)
  trian = get_triangulation(s)
  f = CellData.GenericCellField(object,trian,PhysicalDomain())
  s(f)
end

function Gridap.FESpaces._cell_vals(fs::SingleFieldFESpace,cf_ambient::CellField)
  println("my cell vals")
  s = get_fe_dof_basis(fs)
  s(cf_ambient)
end

function pinv(J::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  @check D1 !== D2
  Jt = transpose(J)
  inv(Jt⋅J)⋅Jt
end

################################################################################
##### tangent of analytical vector field in ambient space
################################################################################
function tangent_f(vec_func::Function)
  function _tangent_f(X)
    vec = vec_func(X)

    normal_vec = 1/sqrt(X[1]*X[1] + X[2]*X[2] + X[3]*X[3])*VectorValue(X[1],X[2],X[3])
    normal_comp = (vec⋅normal_vec)*normal_vec
    tangent_comp = vec - normal_comp

    @assert norm(normal_vec) ≈ 1.0 # check length of normal vector = 1
    @assert dot(normal_comp, tangent_comp) <= 9.9e-15 # check normal and tangent components are perpendicular

    tangent_comp

  end
end


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


vec(x) = VectorValue(-x[2],x[1],0.0)


## interpolate analytic function into ambient space
RT_ambient = FESpace(ambient_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
analytic_u_ambient = interpolate(tangent_f(vec),RT_ambient)

## projection analytic function into ambient space
project_cell_field = map(p->GenericField(u_projection(p,tangent_f(vec))),panel_ids)
uh_ambient_project = interpolate(project_cell_field,RT_ambient)


## interpolate mapped analytic function into parametric space
RT = FESpace(manifold_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
cell_field = map(p->GenericField(u_parametric(p,tangent_f(vec))),panel_ids)
uh = interpolate(cell_field,RT)

## map parametric FEFunction back to ambient space
mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)
inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)

Jt = lazy_map(Broadcasting(gradient),mapping)
J = lazy_map(Operation(transpose),Jt)
pinvJ = lazy_map(Operation(pinv),J)

_uh = change_domain(uh.cell_field,ReferenceDomain(),PhysicalDomain())

_cf_mapped = lazy_map(Broadcasting(⋅),J,get_data(_uh))
cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,inv_mapping)

cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() ) # ambient cell field

### Convert ambient cellfield into a FEFunction by evaluating the field at the dofs
uh_ambient = interpolate(cf_ambient,RT_ambient)

## compute J pinv(J) uh_ambient
change = lazy_map(Broadcasting(⋅),J,pinvJ)
_change = lazy_map(Broadcasting(∘),change,inv_mapping)
cf_change = CellData.GenericCellField(_change,Ω_ambient,PhysicalDomain() )
cf_out = cf_change ⋅cf_ambient
uh_out = interpolate(cf_out,RT_ambient)


writevtk(Ω_ambient,dir*"/ambient_tangent",
cellfields=["f"=>tangent_f(vec),"u_ambient_cf"=>cf_ambient,"uh_ambient"=>uh_ambient,
"eh"=>uh_ambient_project-uh_out,
"u_analytic_ambient"=>analytic_u_ambient,
"projection"=>uh_ambient_project,
"out"=>uh_out ],append=false)


dΩ = Measure(Ω_ambient,2)
l2(uh_ambient_project-uh_out,dΩ)

l2(analytic_u_ambient-uh_ambient,dΩ)
