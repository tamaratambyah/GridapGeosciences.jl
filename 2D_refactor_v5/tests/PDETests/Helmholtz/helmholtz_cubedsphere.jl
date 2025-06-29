using DrWatson
using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots
include("../../../src/initialise.jl")

function uθϕ_scalar(θϕ)
  θ,ϕ = θϕ
  # cos(ϕ)*cos(θ)
  sin(ϕ)
end

# function uX_scalar(X)
#   X[1]*X[2]*X[3]
# end

function uαβ_scalar(p)
  function _h(αβ)
    if p == 2 || p == 3 || p == 4
      return -αβ[1]*αβ[2]
    else
      return αβ[1]*αβ[2]
    end
  end
end

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)


manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

p_fe = 2
degree = 4*(p_fe+1)

################################################################################
##### ambient space
################################################################################
ambient_model = get_ambient_model(manifold_model)
Ω_amb = Triangulation(ambient_model)
dΩ_amb = Measure(Ω_amb,degree)

# cell_field = map(p->GenericField(u_scalar_parametric2ambient(p,uαβ_scalar(p))),panel_ids)
# cell_field = map(p->GenericField(uX_scalar),panel_ids)
 cell_field = map(p->GenericField(u_scalar_latlon2ambient(p,uθϕ_scalar)),panel_ids)
ucf_ambient = CellData.GenericCellField(cell_field,Ω_amb,PhysicalDomain())


# check compatibility in ambient space
# rhs_amb = ucf_ambient + laplacian(ucf_ambient)
# rhs_amb = ucf_ambient
# compat = sum( ∫(rhs_amb)dΩ_amb   )

writevtk(Ω_amb,dir*"/poisson_ambient",cellfields=["u"=>ucf_ambient],append=false)

################################################################################
##### parametric space
################################################################################
####### check compatibility in parametric space
Ω = Triangulation(manifold_model)
m = Metric(cubedsphere,Ω)

dΩ = Measure(Ω,degree)
dΩg =  Measure(m,Ω,degree)

# cell_field = map(p->GenericField(uαβ_scalar(p)),panel_ids)
# cell_field = map(p->GenericField(u_scalar_ambient2parametric(p,uX_scalar)),panel_ids)
cell_field = map(p->GenericField(u_scalar_latlon2parametric(p,uθϕ_scalar)),panel_ids)
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

rhs_amb = ucf_ambient + parametric_cf_2_ambient(manifold_model,surface_laplacian(ucf,m))
rhs_amb(get_cell_points(Ω_amb))./1

mapping = map(x-> PanelRotationField(r1p_3D[x])  ∘ SigmaField(RADIUS) ∘ GnomonicField() , panel_ids)

cf_mapped = lazy_map(Broadcasting(∘),get_data(rhs_amb),mapping)
_rhs = CellData.GenericCellField(cf_mapped,Ω,PhysicalDomain() )
_rhs(get_cell_points(Ω))./1


rhs = ucf + 1.0*(surface_laplacian(ucf,m))
rhs(get_cell_points(Ω))./1
# rhs = ucf
compat = sum( ∫(rhs)dΩ   )
println("Compatibility: ", compat)

_Jcf = map(p->GenericField(Jp(p)),panel_ids)
Jcf = CellData.GenericCellField(_Jcf,Ω,PhysicalDomain())
Gcf = (Operation(transpose)(Jcf)) ⋅ Jcf

# pts = get_cell_points(Ω)

# Jcf(pts)./1

# Gcf(pts)
# m.metric(pts)
# (Gcf - m.metric)(pts)


V = TestFESpace(manifold_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
U = TrialFESpace(V)

stiffnes(u,v) = ∫( (Jcf ⋅ m.inv_metric ⋅ gradient(v))⋅ (Jcf ⋅ m.inv_metric ⋅gradient(u) ) *(m.sq_meas))dΩ
K = assemble_matrix(stiffnes,U,V)
using LinearAlgebra
issymmetric(Array(K))
# map(x->issymmetric.(x),m.inv_metric(get_cell_points(Ω)))


helmholtz_biform(u,v) = ∫( (u*v)* m.sq_meas )dΩ -  stiffnes(u,v)
helmholtz_liform(v) = ∫(  (rhs*v) * m.sq_meas )dΩ
op = AffineFEOperator(helmholtz_biform,helmholtz_liform,U,V)

uh = solve(LUSolver(),op)

# for pid in collect(1:6)
#   mask = (panel_ids.== pid)
#   Ωp = Triangulation(manifold_model,mask)
#   uhp = change_domain(uh,Ωp,ReferenceDomain())
#   m = Metric(cubedsphere,Ωp)

#   ucfp = change_domain(ucf,Ωp,DomainStyle(ucf))
#   writevtk(Ωp,dir*"/poisson_uhp$(pid)",cellfields=["uh"=>uhp,"graduh"=>surface_gradient(uhp,m),
#                                                   "u"=>ucfp,"gradu"=>gradient(ucfp)],append=false)
# end
# using LinearAlgebra
# A = get_matrix(op)
# b = get_vector(op)
# eigvals(Array(A))
# A*ones(size(b))

## map parametric FEFunction back to ambient space
cf_ambient = parametric_cf_2_ambient(manifold_model,uh.cell_field)


#### Compute errors
e = l2(uh-ucf,dΩ)

eg = l2(uh-ucf,dΩg)

e_amb = l2(cf_ambient-ucf_ambient,dΩ_amb)

# rr = parametric_cf_2_ambient(manifold_model,rhs)
writevtk(Ω_amb,dir*"/poisson_ambient",cellfields=["u"=>ucf_ambient,"uh"=>cf_ambient,"e"=>ucf_ambient-cf_ambient],append=false)


errs_p = []
####### solve / plot in individual panels
for pid in collect(1:6)
  mask = (panel_ids.== pid)
  Ω = Triangulation(manifold_model,mask)
  m = Metric(cubedsphere,Ω)

  dΩ = Measure(Ω,degree)
  dΩg =  Measure(m,Ω,degree)

  # cell_field = map(p->GenericField(u_scalar_latlon2parametric(p,uθϕ_scalar)),panel_ids[findall(mask)])
  cell_field = map(p->GenericField(u_scalar_ambient2parametric(p,uX_scalar)),panel_ids[findall(mask)])
  cf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

  rhs = cf + 1.0*(surface_laplacian(cf,m))
  compat = sum( ∫(rhs)dΩ   )
  println("Compatibility: ", compat)
  println( sum(∫( cf)dΩ), "; ", sum(∫( (surface_laplacian(cf,m)))dΩ))

  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U = TrialFESpace(V)

  poisson_biform(u,v) = ∫(u*v)dΩg -  ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
  poisson_liform(v) = ∫(  rhs*v )dΩg
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

  uh = solve(LUSolver(),op)

  #### Compute errors
  e = l2(uh-cf,dΩ)
  eg = l2(uh-cf,dΩg)
  println("Errors: ", e, "; ", eg)
  push!(errs_p,e)
  writevtk(Ω,dir*"/poisson_panel$(pid)",cellfields=["u"=>cf, "uh"=>uh,"e"=>cf-uh],append=false)
end

include("../pde_helpers.jl")

plot()
plot_error(collect(1:6),errs_p)
plot!(framestyle=:box,
xlabel="panel",ylabel="L2(u - uh)")
savefig(plotsdir()*"/helmholtz_cubedsphere")
