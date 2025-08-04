using DrWatson
using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots
include("../../../src/initialise.jl")



function uX_scalar(X)
  X[1]*X[2]*X[3]
end


manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)

manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

p_fe = 2
degree = 4*(p_fe+1)


################################################################################
##### parametric space
################################################################################
####### check compatibility in parametric space
Ω = Triangulation(manifold_model)
m = Metric(cubedsphere,Ω)

dΩ = Measure(Ω,degree)
dΩg =  Measure(m,Ω,degree)

cell_field = map(p->GenericField(u_scalar_ambient2parametric(p,uX_scalar)),panel_ids)
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

pts = get_cell_points(Ω)
ucf(pts)[panel_ids.==1]


lap_ucf = surface_laplacian(ucf,m)

rhs =  lap_ucf
rhs(get_cell_points(Ω))./1
compat = sum( ∫(rhs)dΩ   )
println("Compatibility: ", compat)


V = TestFESpace(manifold_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1, constraint=:zeromean)
U = TrialFESpace(V)

_biform(u,v) = ∫( u*v )dΩg
_liform(v) = ∫( rhs*v )dΩg
op = AffineFEOperator(_biform,_liform,U,V)
rhsh = solve(LUSolver(),op)
# rhsh = interpolate(rhs,U)

#### Compute errors
e = l2(rhs-rhsh,dΩ)

eg = l2(rhs-rhsh,dΩg)

# rr = parametric_cf_2_ambient(manifold_model,rhs)
writevtk(Ω_amb,dir*"/poisson_ambient",cellfields=["u"=>ucf_ambient,"uh"=>cf_ambient,"e"=>ucf_ambient-cf_ambient],append=false)
