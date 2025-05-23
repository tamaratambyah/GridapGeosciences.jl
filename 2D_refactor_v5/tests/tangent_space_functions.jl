using Gridap
include("../src/initialise.jl")
include("../src/Visualization/VisualizationData.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

ambient_model = get_ambient_model(manifold_model)
latlon_model = get_latlon_model(manifold_model)


Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)
Ω_latlon = Triangulation(latlon_model)

pts_parametric = get_cell_points(Ω_parametric)
pts_ambient = get_cell_points(Ω_ambient)
pts_latlon = get_cell_points(Ω_latlon)


################################################################################
##### analytical function in ambient space
################################################################################

function tangent_f(X)
  # VectorValue(-X[2],X[1],0)

  vec = VectorValue(X[1]^2,X[1]*X[2],X[3])
  normal_vec = 1/r*VectorValue(X[1],X[2],X[3])
  normal_comp = (vec⋅normal_vec)*normal_vec

  vec - normal_comp
end

tangent_f(Point(1,1,1))

writevtk(Ω_ambient,dir*"/ambient_tangent",cellfields=["u"=>tangent_f],append=false)


################################################################################
##### mapping:      α,β   -> X,Y,Z ; 2D -> 3D ; J = 3 x 2     ; J^T = 2 x 3
##### mapping_inv:  X,Y,Z -> α,β   ; 3D -> 2D ; J_inv = 2 x 3 ; (J_inv)^T = 3 x 2
##### check (J^T)^{-1} = (J_inv)^T (extra transpose)
##### see that this condition is true
##### see that the determinate of mapping_inv is inf
################################################################################

function Mymeas(Jt::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  J = transpose(Jt)
  d = (det(Jt⋅J))
  if abs(d) < 9.9e-15
    return abs(d)
  else
    return sqrt(d)
  end

end

function MyInvmeas(Jt::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  J = transpose(Jt)
  d = (det(Jt⋅J))
  if abs(d) < 9.9e-15
    return 1/abs(d)
  else
    return 1/sqrt(d)
  end

end


panel1_αβ = (get_grid(manifold_model).parametric_cell_coords)[panel_ids .==1]

for p in collect(1:6)
  cmap = lazy_map(x-> PanelRotationField(r1p_3D[p]) ∘ SigmaField(r) ∘ GnomonicField(), 1:length(panel1_αβ))
  inv_cmap = lazy_map(x-> InvGnomonicField() ∘ InvSigmaField(r)  ∘ PanelRotationField(rp1_3D[p]), 1:length(panel1_αβ))

  latlon_panel1 = lazy_map(GnomonicMap(), panel1_αβ)
  sphere_panel1 = lazy_map(Sigma(),latlon_panel1)
  sphere_X = lazy_map(PanelRotationMap(r1p_3D[p]), sphere_panel1)


  Jt = lazy_map(Broadcasting(∇),cmap)
  Jt_inv = lazy_map(Broadcasting(∇),inv_cmap)

  pinv_Jt = lazy_map(Broadcasting(pinvJt),Jt)

  det_Jt = lazy_map(Operation(meas),Jt)
  det_Jt_inv = lazy_map(Operation(Mymeas),Jt_inv)

  pinv_Jt_x = lazy_map(evaluate,pinv_Jt,panel1_αβ)
  Jt_inv_x = lazy_map(evaluate,Jt_inv,sphere_X)

  det_Jt_x = lazy_map(evaluate,det_Jt,panel1_αβ)
  det_Jt_inv_x = lazy_map(evaluate,det_Jt_inv,sphere_X) ### has infs

  println("p = $p: inverses equal? ", sum(pinv_Jt_x .≈ Jt_inv_x) == length(panel1_αβ),
          "; dets equal? ", sum(det_Jt_x .≈ det_Jt_inv_x) == length(panel1_αβ), )

end

################################################################################
##### Map analytical vector valued function from X,Y,Z -> α,β
##### 1. Map a point α,β -> X,Y,Z
##### 2. Compute u(X,Y,Z)
##### 2. Compute J cooresponding to X,Y,Z -> α,β
##### 2. Evaluate v(α,β) = J^{-T} ⋅ u(X,Y,Z)
################################################################################


function u_ambient_vector(p::Int,uX::Function)
  function _u(αβ)
    latlon_panel1 = evaluate(GnomonicMap(), αβ)
    sphere_panel1 = evaluate(Sigma(),latlon_panel1)
    sphere_panelp = evaluate(PanelRotationMap(r1p_3D[p]), sphere_panel1)

    cmap = PanelRotationField(r1p_3D[p]) ∘ SigmaField(r) ∘ GnomonicField()
    inv_cmap = InvGnomonicField() ∘ InvSigmaField(r)  ∘ PanelRotationField(rp1_3D[p])

    Jt = ∇(cmap)
    Jt_inv = ∇(inv_cmap)

    Jt_x = evaluate(Jt,αβ)
    pinv_Jt_x = evaluate(pinvJt(Jt),αβ)

    det_Jt_x = evaluate(Mymeas,Jt_x)

    (det_Jt_x) * transpose(pinv_Jt_x) ⋅ uX(sphere_panelp)
  end
end






################################################################################
##### Vector valued functions
################################################################################

dΩ = Measure(Ω_parametric,2)
RT = FESpace(manifold_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)

## interpolate everywhere
free_values = zero_free_values(RT)
dirichlet_values = zero_dirichlet_values(RT)

## get cell vals
s = get_fe_dof_basis(RT)
trian = get_triangulation(s)
cell_field = map(p->GenericField(u_ambient_vector(p,tangent_f)),panel_ids)
f = CellData.GenericCellField(cell_field,trian,PhysicalDomain())
cell_vals = s(f)

## interpolate!
gather_free_and_dirichlet_values!(free_values,dirichlet_values,RT,cell_vals)
uh = FEFunction(RT,free_values,dirichlet_values)




# writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>uh],append=false)

mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)
inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)

Jt = lazy_map(Broadcasting(gradient),mapping)
pinv_Jt = lazy_map(Operation(pinvJt),Jt)
det_Jt = lazy_map(Operation(Mymeas),Jt)
change = lazy_map(*,det_Jt,pinv_Jt)

_uh = change_domain(uh,Ω_parametric,PhysicalDomain())

_cf_mapped = lazy_map(Broadcasting(⋅),change,get_data(_uh))
cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,inv_mapping)

cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(pts_ambient)

writevtk(Ω_ambient,dir*"/ambient_tangent",cellfields=["f"=>tangent_f,"u"=>cf_ambient],append=false)
