"""
Manufacture solutions for Possion problem on the cubed sphere mesh.

Solve Δu = -f on Ω, where
  - Δ is the Laplace-Beltrami operator,
  - Ω is the parametric space of manifold M
  - f = -surf_lap(u_ex) is a manufactured rhs

The weak form in terms of surface operators is:
  ∫ (surf_grad(u)) ⋅ (grad(v)) dΩg = ∫ v(-surf_lap(u_ex)) dΩg
where
  - ∫ is in the parametric space
  - surf_grad, surf_lap, dΩg account for the metric
  - grad is the standard flat gradient in the parametric space

The parametric space is Ω = [-π/4,π/4]^2 doubly periodic
  - need to enforce zero mean constraint
Consider piecewise u_ex that is periodic, zeromean and in FE space:
  - u_ex = ( -α*(α + π/4 ) if α < 0 ;
              α*(α - π/4 ) is α > 0 )
Find that the error is large
"""


using Gridap
include("../../src/initialise.jl")

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
# function h(p)
  function h(αβ)
    # cos(4*αβ[1])
    if αβ[1] < 0.0
      return -αβ[1]*(αβ[1] + π/4)
    else
      return αβ[1]*(αβ[1] - π/4)
    end
  end
# end

# cellf = map(p->GenericField(h(p)),panel_ids)
# cf_parametric = CellData.GenericCellField(cellf,Ω_parametric,PhysicalDomain())

cf_parametric = CellField(h,Ω_parametric)
rhs = -1.0*surface_laplacian(cf_parametric,m)


############# Analytic function on ambient space
# hambient(X) = X[1]*X[2]*X[3]

# cellf = map(p->GenericField(u_scalar_ambient2parametric(p,hambient)),panel_ids)
# cf_parametric = CellData.GenericCellField(cellf,Ω_parametric,PhysicalDomain())
# rhs = -1.0*surface_laplacian(cf_parametric,m)


############# FE problem
V = FESpace(manifold_model, ReferenceFE(lagrangian,Float64,p); conformity=:H1,
  constraint=:zeromean)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫( surface_gradient(u,m)⋅gradient(v) )dΩg
poisson_liform(v) =  ∫( (rhs*v) )dΩg

op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(uh-cf_parametric,dΩg)
l2(uh-cf_parametric,dΩ_parametric)

########## Map to ambient space to visualise
ss, uh_ambient_analytic = parametric_cf_2_ambient(manifold_model,degree,cf_parametric)
ss, uh_ambient = parametric_cf_2_ambient(manifold_model,degree,uh)

# plot in ambient space
writevtk(Ω_ambient,dir*"/possion_ambient",
        cellfields=["uh"=>uh_ambient,"u"=>uh_ambient_analytic,
                  "e"=>uh_ambient-uh_ambient_analytic],append=false)
writevtk(Ω_parametric,dir*"/possion_parametric",
        cellfields=["uh"=>uh,"u"=>cf_parametric,
                  "e"=>uh-cf_parametric],append=false)
# writevtk(manifold_model,dir*"/model",append=false)
