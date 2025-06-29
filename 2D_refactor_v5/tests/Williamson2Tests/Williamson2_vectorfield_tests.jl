"""
The purpose of this test is to project a vectorfield from θ,ϕ -> α,β on the
parametric raviart thomas space. Such vectorfields are defined in Williamson1992.
The resulting parametric vector field is mapped back to ambient space and
visualised.
The determinant of the forward and inverse map θ,ϕ -> X,Y,Z is also visualised
to show that the inverse map is ill defined (det J = 0)
"""

using Gridap
include("../../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

ambient_model = get_ambient_model(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)

## map parametric FEFunction back to ambient space
mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(RADIUS) ∘ GnomonicField() , panel_ids)
inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(RADIUS) ∘ PanelRotationField(rp1_3D[x]), panel_ids)

# uθϕ(θϕ) = VectorValue(cos(θϕ[1]),0.0)
# uθϕ(θϕ) = VectorValue(-sin(θϕ[1]),0.0)
# uθϕ(θϕ) = VectorValue(-sin(θϕ[2]),0.0)
# uθϕ(θϕ) = VectorValue(cos(θϕ[1])*sin(θϕ[2]),-sin(θϕ[1]))
uθϕ(θϕ) = VectorValue(cos(θϕ[1])*cos(θϕ[2]),0.0)

### projection machinary
vvec = u_vector_latlon2ambient(uθϕ)

RT_ambient = FESpace(ambient_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
analytic_u_ambient = interpolate(tangent_f(vvec),RT_ambient)

## projection analytic function into ambient space
project_cell_field = map(p->GenericField(u_projection(p,tangent_f(vvec))),panel_ids)
uh_ambient_project = interpolate(project_cell_field,RT_ambient)


## interpolate mapped analytic function into parametric space
RT = FESpace(manifold_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
cell_field = map(p->GenericField(u_vector_ambient2parametric(p,tangent_f(vvec))),panel_ids)
uh = interpolate(cell_field,RT)

## map parametric FEFunction back to ambient space
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


writevtk(Ω_ambient,dir*"/ambient_tangent_latlon",
          cellfields=["f"=>tangent_f(vvec),"u_ambient_cf"=>cf_ambient,"uh_ambient"=>uh_ambient,
          "eh"=>uh_ambient_project-uh_out,
          "u_analytic_ambient"=>analytic_u_ambient,
          "projection"=>uh_ambient_project,
          "out"=>uh_out ],append=false)

dΩ = Measure(Ω_ambient,2)
l2(uh_ambient_project-uh_out,dΩ) ≈ 1e-6
l2(analytic_u_ambient-uh_ambient,dΩ)



#### plot the jacobian of lat lon mapping
function Mymeas(J::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  Jt = transpose(J)
  d = det(Jt⋅J)

  if d > 0
    return sqrt(d)
  else
    return d
  end
end

latlon_mapping = lazy_map(x->SigmaField(RADIUS), panel_ids)
_Jt = lazy_map(Broadcasting(gradient),latlon_mapping)
_J = lazy_map(Operation(transpose),_Jt)
dets = lazy_map(Operation(Mymeas),_J)
dets_mapped = lazy_map(Broadcasting(∘),dets,inv_mapping)
cf_dets = CellData.GenericCellField(dets_mapped,Ω_ambient,PhysicalDomain() )


inv_latlon_mapping = lazy_map(x->InvSigmaField(RADIUS), panel_ids)
inv_Jt = lazy_map(Broadcasting(gradient),inv_latlon_mapping)
inv_J = lazy_map(Operation(transpose),inv_Jt)
inv_dets = lazy_map(Operation(Mymeas),inv_J)
cf_inv_dets = CellData.GenericCellField(inv_dets,Ω_ambient,PhysicalDomain() )

writevtk(Ω_ambient,dir*"/dets",
          cellfields=["det"=>cf_dets, "inv_deet"=>cf_inv_dets ],append=false)
