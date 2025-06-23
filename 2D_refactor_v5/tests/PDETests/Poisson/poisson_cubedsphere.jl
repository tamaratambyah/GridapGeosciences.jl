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

The parametric space is Ω = [-π/4,π/4]^2.
  - need to enforce zero mean constraint
  - need to account for boundary
Consider piecewise u_ex that is periodic, zeromean and in FE space:
  - u_ex = ( -α*(α + π/4 ) if α < 0 ;
              α*(α - π/4 ) is α > 0 )
Find that the error converges
"""

using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots, LaTeXStrings
include("../../../src/initialise.jl")
include("../pde_helpers.jl")


p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

function poisson_cubedsphere_panel(n,p,degree,uex)

  model = CartesianDiscreteModel((-π/4,π/4, -π/4,π/4), (n,n))
  Ω = Triangulation(model)
  m = Metric(cubedsphere,Ω)

  dΩ = Measure(Ω,degree)
  dΩg =  Measure(m,Ω,degree)

  Γ = BoundaryTriangulation(model)
  dΓ = Measure(Γ,degree)
  dΓg = Measure(m,Γ,degree)
  n_Γ = get_normal_vector(Γ)

  ucf = CellField(uex,Ω)

  rhs = -1.0*surface_laplacian(ucf,m)
  h = surface_gradient(ucf,m)⋅n_Γ

  # check compatibility
  compat = sum(∫(ucf)dΩ  ), sum(∫( surface_laplacian(ucf,m))dΩ  )
  println("Compatibility: ", compat)

  #### FE Problem -- zero mean
  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1, constraint=:zeromean)
  U = TrialFESpace(V)

  poisson_biform(u,v) =  ∫( surface_gradient(u,m)⋅gradient(v) )dΩg
  poisson_liform(v) =  ∫( rhs*v )dΩg + ∫( v*h )dΓg
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

  uh = solve(LUSolver(),op)

  #### Compute errors
  e = l2(uh-ucf,dΩ)
  eg = l2(uh-ucf,dΩg)
  return e,eg

end


################################################################################
#### Analytic parametric space for a single panel
################################################################################

global RADIUS = 1.0*sqrt(3.0)


function u1(x)
  if x[2] < 0.0
    return -x[2]*(x[2] + π/4)
  else
    return x[2]*(x[2] - π/4)
  end

end
u2(x) = cos(4*x[1])
u3(x) = cos(4*x[1])*sin(4*x[2])


uex_funcs = Dict{Symbol,Any}()
uex_funcs[:u1] = u1
uex_funcs[:u2] = u2
uex_funcs[:u3] = u3

errs = []
errs_g = []
for (key, val) in uex_funcs
  for n in ns
    e, eg = poisson_cubedsphere_panel(n,p,degree,val)
    push!(errs,e)
    push!(errs_g,eg)
  end
end


leginf = map(x->string(x),collect(keys(uex_funcs)))

#
plot()
plot_error(ns,errs;leginf=leginf,ls=fill(:solid,length(uex_funcs)))
plot_error(ns,errs_g;ls=fill(:dash,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:bottomright)
plot!(ns,dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/poisson_cubedsphere_panel_convergence")
