using DrWatson
using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots
include("../../../src/initialise.jl")

# function uθϕ_scalar(θϕ)
#   θ,ϕ = θϕ
#   # cos(ϕ)*cos(θ)
#   sin(ϕ)
# end

function uX_scalar(X)
  X[1]*X[2]*X[3]
end

function uαβ_scalar(p)
  function _h(αβ)
    α,β = αβ
    return α

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

cell_field = map(p->GenericField(u_scalar_parametric2ambient(p,uαβ_scalar(p))),panel_ids)
# cell_field = map(p->GenericField(uX_scalar),panel_ids)
# cell_field = map(p->GenericField(u_scalar_latlon2ambient(p,uθϕ_scalar)),panel_ids)
ucf_ambient = CellData.GenericCellField(cell_field,Ω_amb,PhysicalDomain())


# check compatibility in ambient space
rhs_amb = ucf_ambient + laplacian(ucf_ambient)
rhs_amb = ucf_ambient
compat = sum( ∫(rhs_amb)dΩ_amb   )

writevtk(Ω_amb,dir*"/poisson_ambient",cellfields=["u"=>ucf_ambient],append=false)

################################################################################
##### parametric space
################################################################################
####### check compatibility in parametric space
Ω = Triangulation(manifold_model)
m = Metric(cubedsphere,Ω)

dΩ = Measure(Ω,degree)
dΩg =  Measure(m,Ω,degree)

cell_field = map(p->GenericField(uαβ_scalar(p)),panel_ids)
# cell_field = map(p->GenericField(u_scalar_ambient2parametric(p,uX_scalar)),panel_ids)
# cell_field = map(p->GenericField(u_scalar_latlon2parametric(p,uθϕ_scalar)),panel_ids)
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

pts = get_cell_points(Ω)
ucf(pts)[panel_ids.==1]


lap_ucf = surface_laplacian(ucf,m)

rhs = ucf + lap_ucf #1.0*(surface_laplacian(ucf,m))
rhs(get_cell_points(Ω))./1
compat = sum( ∫(rhs)dΩ   )
println("Compatibility: ", compat)


V = TestFESpace(manifold_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1,constraint=:zeromean)
U = TrialFESpace(V)

Jp(p) = x -> r1p_3D[p]⋅Jpanel1(x)
_Jcf = map(p->GenericField(Jp(p)),panel_ids)
Jcf = CellData.GenericCellField(_Jcf,Ω,PhysicalDomain())

stiffnes(u,v) = ∫( (Jcf ⋅ m.inv_metric ⋅ gradient(v))⋅ (Jcf ⋅ m.inv_metric ⋅gradient(u) ) *(m.sq_meas))dΩ


helmholtz_biform(u,v) = ∫(  (u*v) * m.sq_meas )dΩ - stiffnes(u,v)
helmholtz_liform(v) = ∫(  (rhs*v) * m.sq_meas )dΩ
op = AffineFEOperator(helmholtz_biform,helmholtz_liform,U,V)

uh = solve(LUSolver(),op)
e = l2(uh-ucf,dΩ)

## map parametric FEFunction back to ambient space
cf_ambient = parametric_cf_2_ambient(manifold_model,uh.cell_field)


#### Compute errors
e = l2(uh-ucf,dΩ)

eg = l2(uh-ucf,dΩg)

e_amb = l2(cf_ambient-ucf_ambient,dΩ_amb)

# rr = parametric_cf_2_ambient(manifold_model,rhs)
writevtk(Ω_amb,dir*"/poisson_ambient",cellfields=["u"=>ucf_ambient,"uh"=>cf_ambient,"e"=>ucf_ambient-cf_ambient],append=false)
