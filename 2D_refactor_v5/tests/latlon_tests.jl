using Gridap
include("../src/initialise.jl")


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
##### Map analytical function from α,β -> θ,ϕ
################################################################################
function u_latlon(p::Int,uθϕ::Function)
  function _u(αβ)
    latlon_panel1 = evaluate(GnomonicMap(), αβ)
    sphere_panel1 = evaluate(Sigma(),latlon_panel1)
    sphere_panelp = evaluate(PanelRotationMap(r1p_3D[p]), sphere_panel1)
    latlon_panelp = evaluate(Sigma(),sphere_panelp)
    uθϕ(latlon_panelp)
  end
end



################################################################################
##### Scalar valued functions
################################################################################
function uθϕ_scalar(θϕ)
  θ,ϕ = θϕ
  ϕ
end

cell_field = map(p->GenericField(u_latlon(p,uθϕ_scalar)),panel_ids)
cf_parametric = CellData.GenericCellField(cell_field,Ω_parametric,PhysicalDomain())

dΩ = Measure(Ω_parametric,2)
H1 = FESpace(manifold_model,ReferenceFE(lagrangian,Float64,1), conformity=:H1)
biform(u,v) = ∫(u*v)dΩ
liform(v) = ∫(v*cf_parametric )dΩ

op = AffineFEOperator(biform,liform,H1,H1)
uh = solve(LUSolver(),op)

writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>uh],append=false)


cdofs = get_cell_dof_values(uh)
for p in collect(1:6)
  panel_vals = cdofs[panel_ids.==p]
  println("---p = $p: Maximum = ", maximum(maximum.(panel_vals)) )
  println("---p = $p: Minimum = ", minimum(minimum.(panel_vals)) )
end

mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)
cf_mapped = lazy_map(Broadcasting(∘),get_data(uh),mapping)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)

ambient_vals = cf_ambient(pts_ambient)

for p in collect(1:6)
  panel_vals = ambient_vals[panel_ids.==p]
  println("---p = $p: Maximum = ", maximum(maximum.(panel_vals)) )
  println("---p = $p: Minimum = ", minimum(minimum.(panel_vals)) )
end


################################################################################
##### Vector valued functions
################################################################################
function uθϕ_vector(θϕ)
  θ,ϕ = θϕ
  VectorValue(0.0,cos(ϕ))
end

dΩ = Measure(Ω_parametric,2)
RT = FESpace(manifold_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)

## interpolate everywhere
free_values = zero_free_values(RT)
dirichlet_values = zero_dirichlet_values(RT)

## get cell vals
s = get_fe_dof_basis(RT)
trian = get_triangulation(s)
cell_field = map(p->GenericField(u_latlon(p,uθϕ_vector)),panel_ids)
f = CellData.GenericCellField(cell_field,trian,PhysicalDomain())
cell_vals = s(f)

## interpolate!
gather_free_and_dirichlet_values!(free_values,dirichlet_values,RT,cell_vals)
uh = FEFunction(RT,free_values,dirichlet_values)

# biform(u,v) = ∫( u⋅v )dΩ
# liform(v) = ∫( f⋅v )dΩ
# op = AffineFEOperator(biform,liform,RT,RT)
# uh = solve(LUSolver(),op)

writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>uh],append=false)


cdofs = get_cell_dof_values(uh)
for p in collect(1:6)
  panel_vals = cdofs[panel_ids.==p]
  println("---p = $p: Maximum = ", maximum(maximum.(panel_vals)) )
  println("---p = $p: Minimum = ", minimum(minimum.(panel_vals)) )
end


mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)
inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)

# _mapping = lazy_map(∘,mapping,cmaps)

# attempt 1
_uh = change_domain(uh.cell_field,Ω_parametric,PhysicalDomain())
_cf_mapped = lazy_map(Broadcasting(push_∇),get_data(_uh),mapping)
cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,inv_mapping)

cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(pts_ambient)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)


### attempt 2
# Inverse of the Jacobian transpose
cell_Jt = lazy_map(∇,mapping)
cell_invJt = lazy_map(Operation(pinvJt),cell_Jt)

xt = pts_parametric.cell_phys_point
cell_Jt_x = lazy_map(evaluate,cell_Jt,xt)
cell_invJt_x = lazy_map(evaluate,cell_invJt,xt)

p2 = cell_invJt_x[panel_ids.==2]
p5 = cell_invJt_x[panel_ids.==5]


_cf_mapped = lazy_map(Broadcasting(⋅),cell_invJt,get_data(uh))
_cvals = lazy_map(evaluate,_cf_mapped,xt)
# cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,inv_mapping)

cf_ambient = CellData.GenericCellField(_cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(pts_ambient)

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)



nf_nodes, nf_moments = ReferenceFEs._RT_nodes_and_moments(Float64,QUAD,1,GenericField(identity))
face_own_dofs = ReferenceFEs._face_own_dofs_from_moments(nf_moments)

fcips, fmoments = ReferenceFEs._RT_face_values(QUAD,Float64,1,GenericField(identity))
