"""
consider the function f(θ,ϕ) = cos(θ), defined in the space of latlon
this is equivalent to f(X,Y,Z) = X/sqrt(X^2 + Y^2), defined on the ambient space of the sphere

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


function uθϕ(θϕ)
  θ,ϕ = θϕ
  cos(θ)
end

function grad_uθϕ(θϕ)
  gradient(uθϕ)(θϕ)
end

function uX(X) # = cos(θ)
  x,y,z = X
  x/sqrt(x^2+y^2)
end

function grad_uX(X)
  gradient(uX)(X)
end

# evaluate analytical function U(θϕ) on parametric space αβ
function u_latlon(p::Int,U::Function) # p = panel id
  function _u(αβ)
    cache = return_cache(GnomonicField(),αβ)
    θϕ1 = evaluate!(cache,GnomonicField(),αβ)
    X1 = evaluate(SigmaField(r),θϕ1)
    Xp = evaluate(PanelRotationField(r1p_3D[p]),X1)
    θϕ = evaluate(SigmaField(r),Xp)

    U(θϕ)
  end
end

# evaluate analytical function U(X,Y,Z) on parametric space αβ
function u_ambient(p::Int,U::Function) # p = panel id
  function _u(αβ)
    cache = return_cache(GnomonicField,αβ)
    θϕ1 = evaluate!(cache,GnomonicField(),αβ)

    cache = return_cache(SigmaField(r),θϕ1)
    X1 = evaluate!(cache,SigmaField(r),θϕ1)

    cache = return_cache(PanelRotationField(r1p_3D[p]),X1)
    Xp = evaluate!(cache,PanelRotationField(r1p_3D[p]),X1)

    U(Xp)
  end
end

## change the domain
cmap_ambient = map(x-> InvGnomonicField() ∘ SigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)  # ambient -> parametric
cmap_parametric = map(x->   PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)  # parametric -> ambient

################################################################################
#### Evaluate analytical functions in ambient space
################################################################################
uX_cf = CellField(uX,Ω_ambient)
uX_cvals = uX_cf(pts_ambient)

grad_uX_cf = CellField(grad_uX,Ω_ambient)
grad_uX_cvals = grad_uX_cf(pts_ambient)

writevtk(Ω_ambient,dir*"/williamson2",
      cellfields=["u"=>uX_cf,"gradu"=>grad_uX_cf],append=false)



################################################################################
#### Velocity function
################################################################################
# compute function in parametric space
cell_field = map(p->GenericField(u_ambient(p,uX)),panel_ids)
cf_parametric = GenericCellField(cell_field,Ω,PhysicalDomain())
cvals = cf_parametric(pts)
sum(cvals .≈ uX_cvals)

# map to ambient
cf_mapped = lazy_map(Broadcasting(∘),get_data(cf),cmap_ambient)
velocity_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = velocity_ambient(get_cell_points(Ω_ambient))

sum(cvals_ambient .≈ uX_cvals) == num_cells(sphere_manifold_model)

################################################################################
#### Gradient function
################################################################################
# compute gradient in parametric space
# grad_cell_field = map(p->GenericField(u_ambient(p,grad_uX)),panel_ids)
# grad_cf = GenericCellField(grad_cell_field,Ω,PhysicalDomain())
grad_cf = gradient(cf_parametric)
cvals_grad = grad_cf(pts)

# map to ambient
cf_mapped = lazy_map(Broadcasting(push_∇),get_data(grad_cf),cmap_ambient)
grad_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_grad_ambient = grad_ambient(pts_ambient)

sum(cvals_grad_ambient .≈ grad_uX_cvals) == num_cells(sphere_manifold_model)



################################################################################
#### Surface metric via the metric
################################################################################

# surface gradient
m = Metric(sphere_manifold_model)
surf_gradcf = surface_gradient(cf_parametric,m)
surf_gradcf(pts)



## map to ambient
cf_mapped = lazy_map(Broadcasting(push_∇),get_data(surf_gradcf),cmap_ambient)
surf_grad_ambient = CellData.similar_cell_field(cf_parametric,cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_surf_grad_ambient = surf_grad_ambient(get_cell_points(Ω_ambient))










### debug

∇map = lazy_map(Broadcasting(∇),cmap_parametric)

# map to ambient
_cf_mapped = lazy_map(Broadcasting(push_∇),get_data(grad_cf),∇map)
cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,cmap_ambient)

grad_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_grad_ambient = grad_ambient(get_cell_points(Ω_ambient))

display(cvals_grad_ambient./1)


_grad = CellData.GenericCellField(_cf_mapped,Ω,PhysicalDomain() )

display(_grad(pts))
