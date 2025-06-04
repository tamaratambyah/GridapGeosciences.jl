using Gridap
using Test
include("../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

ambient_model = get_ambient_model(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)

pts_parametric = get_cell_points(Ω_parametric)
pts_ambient = get_cell_points(Ω_ambient)

## map parametric FEFunction back to ambient space
mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)
inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)


biform(u,v,dΩ) = ∫(v*u)dΩ
liform(v,f,dΩ)   = ∫(v*f)dΩ

################################################################################
##### Map analytical scalar valued function from θ,ϕ -> α,β
##### 1. Map a point α,β -> θ,ϕ
##### 2. Evaluate g(α,β) = f(θ,ϕ)
################################################################################
function u_latlon_scalar(p::Int,uθϕ::Function)
  function _u(αβ)
    latlon_panel1 = evaluate(GnomonicMap(), αβ)
    sphere_panel1 = evaluate(Sigma(),latlon_panel1)
    sphere_panelp = evaluate(PanelRotationMap(r1p_3D[p]), sphere_panel1)
    latlon_panelp = evaluate(Sigma(),sphere_panelp)
    uθϕ(latlon_panelp)
  end
end

function streamfunction(α::Float64)
  function _streamfunction(θϕ)
    θ,ϕ = θϕ
    -(sin(ϕ)*cos(α) - cos(θ)*cos(ϕ)*sin(α))
  end
end

# cell field on parametric space
scalar_cf = map(p->GenericField(u_latlon_scalar(p,streamfunction(0.05))),panel_ids)
cf_parametric = CellData.GenericCellField(scalar_cf,Ω_parametric,PhysicalDomain())

# evaluation at quadrature points
dΩ_parametric = Measure(Ω_parametric,2)
L2_parametric = FESpace(manifold_model,ReferenceFE(lagrangian,Float64,1), conformity=:L2)
op_parametric = AffineFEOperator((u,v)->biform(u,v,dΩ_parametric),
                                  v->liform(v,cf_parametric,dΩ_parametric),
                                  L2_parametric,L2_parametric)
uh_parametric = solve(LUSolver(),op_parametric)

_uh_parametric = change_domain(uh_parametric,ReferenceDomain(),PhysicalDomain())

# map to cell field on ambient space
cf_mapped = lazy_map(Broadcasting(∘),get_data(_uh_parametric),inv_mapping)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cf_ambient(pts_ambient)

# gradient(cf_ambient)(pts_ambient)


# evaluation at quadrature points
dΩ_ambient = Measure(Ω_ambient,2)
L2_ambient = FESpace(ambient_model,ReferenceFE(lagrangian,Float64,1), conformity=:L2)
op_ambient = AffineFEOperator((u,v)->biform(u,v,dΩ_ambient),
                              v->liform(v,cf_ambient,dΩ_ambient),
                              L2_ambient,L2_ambient)
uh_ambient = solve(LUSolver(),op_ambient)



##### surface gradient
m = Metric(manifold_model)
surf_grad = surface_gradient(uh_parametric,m)

_surf_grad = change_domain(surf_grad,ReferenceDomain(),PhysicalDomain())

## map parametric FEFunction to ambient space
Jt = lazy_map(Broadcasting(gradient),mapping)
J = lazy_map(Operation(transpose),Jt)
pinvJ = lazy_map(Operation(pinv),J)

_surf_grad_mapped = lazy_map(Broadcasting(⋅),J,get_data(_surf_grad))
surf_grad_mapped = lazy_map(Broadcasting(∘),_surf_grad_mapped,inv_mapping)

surf_grad_ambient = CellData.GenericCellField(surf_grad_mapped,Ω_ambient,PhysicalDomain() ) # ambient cell field

l2(surf_grad_ambient-gradient(uh_ambient),dΩ_ambient)
# l2(surf_grad_ambient-gradient(cf_ambient),dΩ_ambient)


writevtk(Ω_ambient,dir*"/ambient_streamfunc",
    cellfields=["u"=>cf_ambient,"surf_grad"=>surf_grad_ambient,
                # "grad_ambient"=>gradient(cf_ambient),
                "grad_ambient"=>gradient(uh_ambient),
                "e"=>surf_grad_ambient-gradient(uh_ambient)],append=false)
