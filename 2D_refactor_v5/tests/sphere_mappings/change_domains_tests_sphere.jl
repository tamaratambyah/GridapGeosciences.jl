using Gridap
include("../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
# manifold_model = Adaptivity.refine(manifold_model)

ambient_model = get_ambient_model(manifold_model)
latlon_model = get_latlon_model(manifold_model)

panel_ids = get_panel_ids(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)
Ω_latlon = Triangulation(latlon_model)

pts_parametric = get_cell_points(Ω_parametric)
pts_ambient = get_cell_points(Ω_ambient)
pts_latlon = get_cell_points(Ω_latlon)



cmap_ambient = map(x-> InvGnomonicField() ∘ SigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)  # ambient -> parametric
cmap_parametric = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)  # parametric -> ambient


################################################################################
#### Scalar-valued functions in the parametric space
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

quad_parametric = CellQuadrature(Ω_parametric,2)
quad_pts_parametric = quad_parametric.cell_point
quad_ambient = CellQuadrature(Ω_ambient,2)
quad_pts_ambient = quad_ambient.cell_point

cc = get_cell_map(ambient_model)
pp = lazy_map(∘,cmap_ambient,cc)

Jt = lazy_map((∇),pp)
det_Jt = lazy_map(Operation(meas),(Jt))
det_Jtx = lazy_map(evaluate,det_Jt,quad_pts_ambient)

cvals_parametric = cf_parametric(quad_pts_parametric)
cvals_ambient = cf_ambient(quad_pts_ambient)




#####
quad_parametric = CellQuadrature(Ω_parametric,2)
quad_pts_parametric = quad_parametric.cell_point
quad_ambient = CellQuadrature(Ω_ambient,2)
quad_pts_ambient = quad_ambient.cell_point

H1 = FESpace(manifold_model,ReferenceFE(lagrangian,Float64,1),conformity=:H1)
alpha(x) = x[1]

uh = interpolate_everywhere(alpha,H1)
uhc = get_data(uh)
cvals_parametric = uh(get_cell_points(quad_parametric))

M =  map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)
Minv = map(x-> InvGnomonicField() ∘ SigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)


cf_mapped = lazy_map(Broadcasting(∘),get_data(uh),Minv)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,DomainStyle(uh) )
cvals_ambient = cf_ambient(get_cell_points(quad_ambient))

Jt = lazy_map((∇),gg)
det_Jt = lazy_map(Operation(meas),(Jt))
det_Jtx = lazy_map(evaluate,det_Jt,pts_ambient)
det_Jt[1](Point(0,0))

# cf_parametric = CellField(alpha,Ω_parametric)
# cvals_parametric = cf_parametric(pts_parametric)

# cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_parametric),cmap_ambient)
# cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
# cvals_ambient = cf_ambient(pts_ambient)

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
#### Vector-valued functions in the parametric space
################################################################################
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
cmap_ambient = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]) , panel_ids)  # ambient -> parametric

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

RT = FESpace(manifold_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
# cell_field = map(p->GenericField(u_latlon(p,uθϕ)),panel_ids)
uh = interpolate_everywhere(u_latlon(1,uθϕ), RT)

cmap_parametric = map(x->  PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)  # parametric -> ambient
Jt = lazy_map(Broadcasting(∇),cmap_parametric)
Jt_inv = lazy_map(Operation(pinvJt),Jt)

# do the piola map, then map to ambient space
cf_mapped =  lazy_map(Broadcasting(⋅),Jt_inv,get_data(uh))
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,DomainStyle(uh))

dΩ_ambient  = CellQuadrature(Ω_ambient,2)
quad_pts_ambient = get_cell_points(dΩ_ambient)

cvals_ambient = cf_ambient(quad_pts_ambient)
writevtk(Ω_ambient,dir*"/ambient_vector",cellfields=["u"=>cf_ambient],append=false)







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
#### Vectorfields in the latlon space -- via interpolation
################################################################################

RT = FESpace(latlon_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
function uθϕ(θϕ)
  θ,ϕ = θϕ
  VectorValue(cos(θ),0.0)
end

uh = interpolate(uθϕ,RT)
writevtk(Ω_latlon,dir*"/latlon_vector",cellfields=["u"=>uh],append=false)


cmap = map(x-> SigmaField(r), panel_ids)
invcmap = map(x-> InvSigmaField(r), panel_ids)


Jt = map(x-> ∇(SigmaField(r)), panel_ids)
Jt_inv = lazy_map(Operation(pinvJt),Jt)

# do the piola map, then map to ambient space
cf_mapped =  lazy_map(Broadcasting(⋅),Jt_inv,get_data(uh))
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,DomainStyle(uh))

dΩ_ambient  = CellQuadrature(Ω_ambient,2)
quad_pts_ambient = get_cell_points(dΩ_ambient)

cvals_ambient = cf_ambient(quad_pts_ambient)
writevtk(Ω_ambient,dir*"/ambient_vector",cellfields=["u"=>cf_ambient],append=false)





_RT = FESpace(ambient_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
free_values = zero_free_values(_RT)
s = get_fe_dof_basis(_RT)
cell_vals =  s(cf_mapped)
gather_free_values!(free_values,_RT,cell_vals)





writevtk(Ω_parametric,dir*"/parametric_vector",cellfields=["u"=>cf_parametric],append=false)
writevtk(Ω_ambient,dir*"/ambient_vector",cellfields=["u"=>cf_ambient],append=false)



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


# cell_field = map(p->GenericField(_uθϕ(p)),panel_ids)
# u_latlon = CellData.GenericCellField(cell_field,Ω_latlon,PhysicalDomain())
# u_latlon = CellField(uθϕ,Ω_latlon)

# cmap = map(x-> SigmaField(r), panel_ids)
# do the piola map, then map to ambient space
# _cf_mapped = lazy_map(Broadcasting(push_∇),get_data(u_latlon),cmap)
# cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,cmap)

# cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
# cvals_ambient = cf_ambient(pts_ambient)

# writevtk(Ω_ambient,dir*"/ambient_vector",cellfields=["u"=>cf_ambient],append=false)

#### Testing of sigma Field
panel1 = get_cell_coordinates(latlon_model)[panel_ids .== 1]
panel2 = get_cell_coordinates(latlon_model)[panel_ids .== 2]
panel3 = get_cell_coordinates(latlon_model)[panel_ids .== 3]
panel4 = get_cell_coordinates(latlon_model)[panel_ids .== 4]
panel5 = get_cell_coordinates(latlon_model)[panel_ids .== 5]
panel6 = get_cell_coordinates(latlon_model)[panel_ids .== 6]

lon = lazy_map(ExtractVectorValues(1),panel2)
lat = lazy_map(ExtractVectorValues(2),panel2)

include("../src/Fields/Sigma_v2.jl")

Jt = ∇(SigmaField2(r))
_Jt = lazy_map(Jt,panel6)


pinvx = lazy_map(Broadcasting(pinvJt),_Jt)

_pinvx = pinvx[4]
println(maximum.(_pinvx))
for i in 1:length(_pinvx)
  println(_pinvx[i])
end


X2 =lazy_map(SigmaField1(r),panel2)
_X2 = lazy_map(SigmaField2(r),panel2)

display(X2[2])
display(_X2[2])



X(θϕ) = cos(θϕ[1])*cos(θϕ[2])
Y(θϕ) = sin(θϕ[1])*cos(θϕ[2])
Z(θϕ) = sin(θϕ[2])


_X(θϕ) = -cos(θϕ[1])*cos( θϕ[2] )
_Y(θϕ) = -sin( θϕ[2] )
_Z(θϕ) = -sin(θϕ[1])*cos( θϕ[2])


n = 5
lats = collect(range(-π/2,π/2,n))[randperm(n)]
lons = collect(range(0.0,2*π,n))[randperm(n)]

θϕs = []
for i in 1:n
  push!(θϕs,Point(lons[i],lats[i]))
end

for i in 1:n
  println(i)
  println(X(θϕs[i]), " ",Y(θϕs[i])," ",Z(θϕs[i]))
  println(_X(θϕs[i])," ",_Y(θϕs[i])," ",_Z(θϕs[i]))
end


X2 =lazy_map(SigmaField1(r),θϕs)
_X2 = lazy_map(SigmaField2(r),θϕs)

display(X2[2])
display(_X2[2])
