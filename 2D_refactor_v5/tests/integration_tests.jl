"""
consider the function f(α,β) = α, defined in the parametric space

want to compute the gradient of this function.
"""

using Gridap
include("../src/initialise.jl")

sphere_manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
# sphere_manifold_model = Adaptivity.refine(sphere_manifold_model)


r = sqrt(3.0)

panel_ids = get_panel_ids(sphere_manifold_model)
Ω = Triangulation(sphere_manifold_model) # parametric space α,β
pts = get_cell_points(Ω)

ambient_model = get_ambient_model(sphere_manifold_model) # ambient sphere X,Y,Z
Ω_ambient = Triangulation(ambient_model)
pts_ambient = get_cell_points(Ω_ambient)


function uαβ(αβ)
  α,β = αβ
  α
end

function grad_uαβ(αβ)
  gradient(uαβ)(αβ)
end

# evaluate analytical function u on parametric space αβ from ambient point
function u_ambient(p::Int,uαβ::Function) # p = panel id
  function _u(Xp)
    cache = return_cache(PanelRotationField(rp1_3D[p]),Xp)
    X1 = evaluate!(cache,PanelRotationField(rp1_3D[p]),Xp)

    cache = return_cache(SigmaField(RADIUS),X1)
    θϕ1 = evaluate!(cache,SigmaField(RADIUS),X1)


    cache = return_cache(InvGnomonicField,θϕ1)
    αβ = evaluate!(cache,GnomonicField(),θϕ1)

    uαβ(αβ)
  end
end

################################################################################
#### Evaluate analytical functions in parametric space
################################################################################
u_cf = CellField(uαβ,Ω)
u_cvals = u_cf(pts)

grad_u_cf = CellField(grad_uαβ,Ω)
grad_u_cvals = grad_u_cf(pts)

writevtk(Ω,dir*"/williamson2",
      cellfields=["u"=>u_cf,"gradu"=>grad_u_cf],append=false)

gradient(u_cf)(pts)

################################################################################
#### Evaluate gradient in ambient space
################################################################################
cmap_ambient = map(x-> InvGnomonicField() ∘ SigmaField(RADIUS) ∘ PanelRotationField(rp1_3D[x]), panel_ids)  # ambient -> parametric
cmap_parametric = map(x->   PanelRotationField(r1p_3D[x]) ∘ SigmaField(RADIUS) ∘ GnomonicField() , panel_ids)  # parametric -> ambient

cf_mapped = lazy_map(Broadcasting(push_∇),get_data(grad_u_cf),cmap_parametric)
cf_ambient = CellData.similar_cell_field(grad_u_cf,cf_mapped ,get_triangulation(grad_u_cf),DomainStyle(grad_u_cf))
cvals_ambient = cf_ambient(get_cell_points(Ω))

################################################################################
#### Surface gradient via the metric
################################################################################

# surface gradient
m = Metric(sphere_manifold_model)
surf_gradcf = surface_gradient(u_cf,m)
surf_gradcf(pts)


## change the domain
# cmap = map(x-> InvGnomonicField() ∘ SigmaField(RADIUS) ∘ PanelRotationField(rp1_3D[x]), panel_ids)  # ambient -> parametric
cmap = map(x->   PanelRotationField(r1p_3D[x]) ∘ SigmaField(RADIUS) ∘ GnomonicField() , panel_ids)  # parametric -> ambient

g = lazy_map(Broadcasting(push_∇),get_data(surf_gradcf),cmap)
gradu = CellData.similar_cell_field(u_cf,g,get_triangulation(u_cf),DomainStyle(u_cf))
display(gradu(pts))
