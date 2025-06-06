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

############# Surface metric and quadrature
p = 2
degree = 2*(p+1)
m = Metric(cubedsphere,Ω_parametric)
dΩg = Measure(m,Ω_parametric,degree)
dΩ_parametric = Measure(Ω_parametric,degree)

############# Analytic function on parametric space
function h(p)
  function _h(αβ)
    if p == 2 || p == 3 || p == 4
      return -αβ[1]*αβ[2]
    else
      return αβ[1]*αβ[2]
    end
  end
end

cellf = map(p->GenericField(h(p)),panel_ids)
cf_parametric = CellData.GenericCellField(cellf,Ω_parametric,PhysicalDomain())
rhs = -1.0*surface_laplacian(cf_parametric,m)


############# Analytic function on ambient space
hambient(X) = X[1]*X[2]*X[3]

cellf = map(p->GenericField(u_scalar_ambient2parametric(p,hambient)),panel_ids)
cf_parametric = CellData.GenericCellField(cellf,Ω_parametric,PhysicalDomain())
rhs = -1.0*surface_laplacian(cf_parametric,m)



############# FE problem


V = FESpace(manifold_model, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫(  surface_gradient(v,m) ⋅ surface_gradient(u,m)  )dΩg
poisson_liform(v)   = ∫( v*rhs )dΩg

op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(uh-cf_parametric,dΩg)
l2(uh-cf_parametric,dΩ_parametric)

########## Map to ambient space to visualise
ss, uh_ambient_analytic = parametric_cf_2_ambient(manifold_model,degree,cf_parametric)
ss, uh_ambient = parametric_cf_2_ambient(manifold_model,degree,uh)

# plot in parametric space
writevtk(Ω_ambient,dir*"/possion_ambient",
        cellfields=["uh"=>uh_ambient,"u"=>uh_ambient_analytic,
                  "e"=>uh_ambient-uh_ambient_analytic],append=false)
