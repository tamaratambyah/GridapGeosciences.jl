"""
The purpose of this test is to map a scalar field from θ,ϕ -> α,β on the
parametric L2 space. Such scalar fields are defined in Williamson1992.
The resulting parametric scalar field is mapped back to ambient space and
visualised.
The minimum and maximum vals on compared panel-wise between uh_parametric and
cf_ambient
"""

using Gridap
include("../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

ambient_model = get_ambient_model(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)

pts_parametric = get_cell_points(Ω_parametric)
pts_ambient = get_cell_points(Ω_ambient)


function uθϕ_scalar(θϕ)
  θ,ϕ = θϕ
  # rem2pi(θ,RoundNearest)
  cos(θ)
end

cell_field = map(p->GenericField(u_scalar_latlon2parametric(p,uθϕ_scalar)),panel_ids)
cf_parametric = CellData.GenericCellField(cell_field,Ω_parametric,PhysicalDomain())

uh_parametric, uh_ambient = parametric_cf_2_ambient(manifold_model,2,cf_parametric)

### compare values
cdofs = get_cell_dof_values(uh_parametric)
ambient_vals = get_cell_dof_values(uh_ambient)
for p in collect(1:6)
  parametric_panel_vals = cdofs[panel_ids.==p]
  ambient_panel_vals = ambient_vals[panel_ids.==p]

  @test maximum(maximum.(parametric_panel_vals)) ≈ maximum(maximum.(ambient_panel_vals))
  @test minimum(minimum.(parametric_panel_vals)) ≈ minimum(minimum.(ambient_panel_vals))
end



################################################################################
#### low level test
################################################################################
## map parametric FEFunction back to ambient space
mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(RADIUS) ∘ GnomonicField() , panel_ids)
inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(RADIUS) ∘ PanelRotationField(rp1_3D[x]), panel_ids)


dΩ = Measure(Ω_parametric,2)
L2 = FESpace(manifold_model,ReferenceFE(lagrangian,Float64,1), conformity=:L2)
biform(u,v) = ∫(u*v)dΩ
liform(v) = ∫(v*cf_parametric )dΩ

op = AffineFEOperator(biform,liform,L2,L2)
uh = solve(LUSolver(),op)

writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>uh],append=false)

## map parametric FEFunction back to ambient space
_uh = change_domain(uh,ReferenceDomain(),PhysicalDomain())

cf_mapped = lazy_map(Broadcasting(∘),get_data(_uh),inv_mapping)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )

writevtk(cubedsphere,Ω_ambient,dir*"/ambient_quad",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_ambient,dir*"/ambient_scalar",cellfields=["u"=>cf_ambient],append=false)

### compare values
cdofs = get_cell_dof_values(uh)
ambient_vals = cf_ambient(pts_ambient)
for p in collect(1:6)
  parametric_panel_vals = cdofs[panel_ids.==p]
  ambient_panel_vals = ambient_vals[panel_ids.==p]

  @test maximum(maximum.(parametric_panel_vals)) ≈ maximum(maximum.(ambient_panel_vals))
  @test minimum(minimum.(parametric_panel_vals)) ≈ minimum(minimum.(ambient_panel_vals))

  println("Parametric ---p = $p: Maximum = ", maximum(maximum.(parametric_panel_vals)),
                                "; Minimum = ", minimum(minimum.(parametric_panel_vals)) )

  println("Ambient ------p = $p: Maximum = ", maximum(maximum.(ambient_panel_vals)),
                            "; Minimum = ", minimum(minimum.(ambient_panel_vals)) )
end
