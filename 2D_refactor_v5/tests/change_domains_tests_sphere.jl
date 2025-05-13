using Gridap
include("../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)

ambient_model = get_ambient_model(manifold_model)
latlon_model = get_latlon_model(manifold_model)

panel_ids = get_panel_ids(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)
Ω_latlon = Triangulation(latlon_model)

pts_parametric = get_cell_points(Ω_parametric)
pts_ambient = get_cell_points(Ω_ambient)
pts_latlon = get_cell_points(Ω_latlon)




################################################################################
#### functions in the parametric space
################################################################################
cmap_ambient = map(x-> InvGnomonicField() ∘ SigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)  # ambient -> parametric


alpha(x) = x[1]

cf_parametric = CellField(alpha,Ω_parametric)
cvals_parametric = cf_parametric(pts_parametric)

cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_parametric),cmap_ambient)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(pts_ambient)

display(cvals_parametric)
display(cvals_ambient)
cvals_parametric ≈ cvals_ambient

#####
beta(x) = x[2]

cf_parametric = CellField(beta,Ω_parametric)
cvals_parametric = cf_parametric(pts_parametric)

cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_parametric),cmap_ambient)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(pts_ambient)

display(cvals_parametric)
display(cvals_ambient)
cvals_parametric ≈ cvals_ambient



writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)



################################################################################
#### Vectorfields in the parametric space
################################################################################
cmap_ambient = map(x-> InvGnomonicField() ∘ SigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)  # ambient -> parametric
cmap_parametric = map(x->   PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)  # parametric -> ambient

function u(αβ)
  α,β = αβ
  VectorValue((α+π/4)*(α-π/4),(β+π/4)*(β-π/4))
  # VectorValue(2*α,β)
end
cf_parametric = CellField(u,Ω_parametric)
cvals_parametric = cf_parametric(pts_parametric)

# do the piola map, then map to ambient space
_cf_mapped = lazy_map(Broadcasting(push_∇),get_data(cf_parametric),cmap_parametric)
cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,cmap_ambient)

cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(pts_ambient)


writevtk(Ω_parametric,dir*"/parametric_vector",cellfields=["u"=>cf_parametric],append=false)
writevtk(Ω_ambient,dir*"/ambient_vector",cellfields=["u"=>cf_ambient],append=false)

display(cvals_parametric)
display(cvals_ambient)


################################################################################
#### Vectorfields in the ambient space
################################################################################
function uX(X)
  x,y,z = X
  Theta = atan(y,x)
  _Theta = rem2pi(atan(y,x),RoundToZero)
  VectorValue(cos(_Theta),0.0)
end

cf_ambient = CellField(uX,Ω_ambient)
writevtk(Ω_ambient,dir*"/ambient_vector",cellfields=["u"=>cf_ambient],append=false)
cf_ambient(pts_ambient)



################################################################################
#### Vectorfields in the latlon space
################################################################################
cmap_ambient = map(x-> InvGnomonicField() ∘ SigmaField(r) ∘ PanelRotationField(rp1_3D[x]) , panel_ids)  # ambient -> parametric
cmap_parametric = map(x->  SigmaField(r) ∘  PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)  # parametric -> ambient

function uθϕ(θϕ)
  θ,ϕ = θϕ
  VectorValue(cos(θ),0.0)
end

function u_latlon(p::Int,uθϕ::Function)
  function _u(αβ)
    latlon_panel1 = evaluate(GnomonicMap(), αβ)
    sphere_panel1 = evaluate(Sigma(),latlon_panel1)
    sphere_panelp = evaluate( PanelRotationMap(rp1_3D[p]), sphere_panel1)
    latlon_panelp = evaluate(Sigma(),sphere_panelp)
    uθϕ(latlon_panelp)
  end
end

cell_field = map(p->GenericField(u_latlon(p,uθϕ)),panel_ids)
cf_parametric = CellData.GenericCellField(cell_field,Ω_parametric,PhysicalDomain())
cvals_parametric = cf_parametric(pts_parametric)

# do the piola map, then map to ambient space
_cf_mapped = lazy_map(Broadcasting(push_∇),get_data(cf_parametric),cmap_parametric)

_cf = CellData.GenericCellField(_cf_mapped,Ω_parametric,PhysicalDomain())
_cf(pts_parametric)


cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,cmap_ambient)

cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(pts_ambient)


writevtk(Ω_parametric,dir*"/parametric_vector",cellfields=["u"=>cf_parametric],append=false)
writevtk(Ω_ambient,dir*"/ambient_vector",cellfields=["u"=>cf_ambient],append=false)

display(cvals_parametric)
display(cvals_ambient)



################################################################################
#### Vectorfields in the latlon space
################################################################################


uθϕ(θϕ) = VectorValue(cos(θϕ[1]),0.0)
function _uθϕ(p)
  function u(θϕ)
    # if p == 2 || p == 5
    #   println("panel 2")
    #  return VectorValue(0.0,0.0)
    # else
    #   θ,ϕ = θϕ
    #   return VectorValue(cos(θ),0.0)
    # end
    uθϕ(θϕ)
  end
end


cell_field = map(p->GenericField(_uθϕ(p)),panel_ids)
u_latlon = CellData.GenericCellField(cell_field,Ω_latlon,PhysicalDomain())
# u_latlon = CellField(uθϕ,Ω_latlon)

cmap = map(x-> SigmaField(r), panel_ids)
# do the piola map, then map to ambient space
_cf_mapped = lazy_map(Broadcasting(push_∇),get_data(u_latlon),cmap)
cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,cmap)

cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(pts_ambient)

writevtk(Ω_ambient,dir*"/ambient_vector",cellfields=["u"=>cf_ambient],append=false)


panel2 = get_cell_coordinates(latlon_model)[panel_ids .== 2]

lon = lazy_map(ExtractVectorValues(1),panel2)
lat = lazy_map(ExtractVectorValues(2),panel2)


Jt = ∇(SigmaField(r))
_Jt = lazy_map(Jt,panel2)

Jtx = _Jt[1][1]
Jx = transpose(Jtx)
transpose(inv(Jtx⋅Jx)⋅Jtx)


pinvx = lazy_map(Broadcasting(pinvJt),_Jt)

_pinvx = pinvx[1]
for i in 1:length(_pinvx)
  println(_pinvx[i])
end
