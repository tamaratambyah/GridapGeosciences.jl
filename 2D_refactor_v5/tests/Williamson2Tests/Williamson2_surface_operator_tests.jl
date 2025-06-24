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

m = Metric(manifold_model)

function streamfunction(α::Float64)
  function _streamfunction(θϕ)
    θ,ϕ = θϕ
    -(sin(ϕ)*cos(α) - cos(θ)*cos(ϕ)*sin(α))
  end
end

# ambient_func(X) = X[1]*X[2]*X[3]

errs = []
for (i,α) in enumerate([0.0,0.05,π/2-0.05,π/2])

  # cell field on parametric space
  scalar_cf = map(p->GenericField(u_scalar_latlon2parametric(p,streamfunction(α))),panel_ids)
  # scalar_cf = map(p->GenericField(u_scalar_ambient2parametric(p,ambient_func)),panel_ids)
  cf_parametric = CellData.GenericCellField(scalar_cf,Ω_parametric,PhysicalDomain())

  uh_parametric, uh_ambient = parametric_cf_2_ambient(manifold_model,2,cf_parametric)


  ##### surface gradient
  surf_grad = surface_gradient(uh_parametric,m)
  surf_grad_ambient = parametric_cf_2_ambient_vector(manifold_model,surf_grad)

  dΩ_ambient = Measure(Ω_ambient,2)
  push!(errs,l2(surf_grad_ambient-gradient(uh_ambient),dΩ_ambient))
  # l2(surf_grad_ambient-gradient(cf_ambient),dΩ_ambient)


  writevtk(Ω_ambient,dir*"/ambient_streamfunc_$i",
      cellfields=["u"=>uh_ambient,"surf_grad"=>surf_grad_ambient,
                  # "grad_ambient"=>gradient(cf_ambient),
                  "grad_ambient"=>gradient(uh_ambient),
                  "e"=>surf_grad_ambient-gradient(uh_ambient)],append=false)
end
println(errs)


################################################################################
#### low level test
################################################################################
α = 0.0
## map parametric FEFunction back to ambient space
mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(RADIUS) ∘ GnomonicField() , panel_ids)
inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(RADIUS) ∘ PanelRotationField(rp1_3D[x]), panel_ids)


scalar_cf = map(p->GenericField(u_scalar_latlon2parametric(p,streamfunction(α))),panel_ids)
cf_parametric = CellData.GenericCellField(scalar_cf,Ω_parametric,PhysicalDomain())


# evaluation at quadrature points
dΩ_parametric = Measure(Ω_parametric,2)
Space_parametric = FESpace(manifold_model,ReferenceFE(lagrangian,Float64,1), conformity=:H1)
op_parametric = AffineFEOperator((u,v)->biform(u,v,dΩ_parametric),
                                  v->liform(v,cf_parametric,dΩ_parametric),
                                  Space_parametric,Space_parametric)
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
surf_grad = surface_gradient(uh_parametric,m)
_surf_grad = change_domain(surf_grad,ReferenceDomain(),PhysicalDomain())

## map parametric FEFunction to ambient space
Jt = lazy_map(Broadcasting(gradient),mapping)
J = lazy_map(Operation(transpose),Jt)
pinvJ = lazy_map(Operation(pinv),J)

_surf_grad_mapped = lazy_map(Broadcasting(⋅),J,get_data(_surf_grad))
surf_grad_mapped = lazy_map(Broadcasting(∘),_surf_grad_mapped,inv_mapping)

surf_grad_ambient = CellData.GenericCellField(surf_grad_mapped,Ω_ambient,PhysicalDomain() ) # ambient cell field


kvec(X) = VectorValue(0.0,0.0,1.0)
ck = CellField(kvec,Ω_ambient)

cross(ck,surf_grad_ambient)(pts_ambient)

### velocity
function  uθϕ(α)
  function uθϕ(θϕ)
    θ,ϕ = θϕ
    u = cos(ϕ)*cos(α) + cos(θ)*sin(ϕ)*sin(α)
    v = -sin(θ)*sin(α)
    VectorValue(u,v)
  end
end
### projection machinary
vvec = u_vector_latlon2ambient(uθϕ(α))

project_cell_field = map(p->GenericField(u_projection(p,tangent_f(vvec))),panel_ids)
vf = CellData.GenericCellField(project_cell_field,Ω_ambient,PhysicalDomain())

l2(surf_grad_ambient-gradient(uh_ambient),dΩ_ambient)
l2(surf_grad_ambient-vvec,dΩ_ambient)
# l2(surf_grad_ambient-gradient(cf_ambient),dΩ_ambient)

function u1(α)
  function _u(θϕ)
    θ,ϕ = θϕ
    cos(ϕ)*cos(α) + sin(ϕ)*cos(θ)*sin(α)
  end
end

scf1 = map(p->GenericField(u_scalar_latlon2ambient(p,u1(α))),panel_ids)
cf1 = CellData.GenericCellField(scf1,Ω_ambient,PhysicalDomain())

function u2(α)
  function _u(θϕ)
    θ,ϕ = θϕ
    -sin(θ)*sin(α)
  end
end

scf2 = map(p->GenericField(u_scalar_latlon2ambient(p,u2(α))),panel_ids)
cf2 = CellData.GenericCellField(scf2,Ω_ambient,PhysicalDomain())

writevtk(cubedsphere,Ω_ambient,dir*"/ambient_streamfunc",
    cellfields=["u"=>uh_ambient,"surf_grad"=>surf_grad_ambient,
                # "grad_ambient"=>gradient(cf_ambient),
                "grad_ambient"=>gradient(uh_ambient),
                "e"=>surf_grad_ambient-gradient(uh_ambient),
                "velocity"=>vf,
                "cross_prod"=>cross(ck,surf_grad_ambient),
                "u1"=>cf1,"u2"=>cf2],append=false)
