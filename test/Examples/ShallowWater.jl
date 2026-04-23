# # Shallow water on the cubed sphere manifold
#
# This example solves the shallow water equations.
# The evolution of the prognostic variables, $\widetilde{\boldsymbol{u}}$ and $\widetilde{\varphi}$ is:
#
# ```math
# \begin{align*}
# \partial_t \widetilde{\boldsymbol{u}}
# + \widetilde{q} \widetilde{\boldsymbol{F}}^\dagger
# + \nabla_{\gamma} \widetilde{\Phi}  &=  0  \quad \text{in} \quad \gamma\times[0,T],
# \\
# \partial_t \widetilde{\varphi}
# + \nabla_{\gamma}\cdot \widetilde{\boldsymbol{F}} &= 0
# \quad  \quad \text{in} \quad \gamma\times[0,T],
# \end{align*}
# ```
#
# The diagnostic variables solve algebraic constraints:
# ```math
# \begin{align*}
# \widetilde{\boldsymbol{F}} &= \widetilde{\varphi}\widetilde{\boldsymbol{u}} , \\
# \widetilde{\Phi} &= \frac{1}{2} \widetilde{\boldsymbol{u}}\cdot\widetilde{\boldsymbol{u}} + g_r \widetilde{\varphi}  ,\\
# \widetilde{q} &= \frac{1 }{\widetilde{\varphi}}(\nabla_{\gamma}^\dagger\cdot \widetilde{\boldsymbol{u}}  + f) ,
# \end{align*}
# ```
#
# where $\gamma$ is the cubed sphere manifold,  $g_r$ is the acceleration due to gravitiy, and $f$ is the Coriolis force.
# We use compatible finite element pairings to solve in the parametric space of the cubed sphere.
# This problem is a differential algebraic equation (DAE), which is solved using a bespoke
# ``DAEFEOperator`` that builds on Gridap's ODE API.

# ## Set up
# First load all required pacakges. In this example, we will use a distributed model
# and iterative solvers. So we initialise MPI.
using GridapGeosciences
using Gridap
using GridapSolvers
using GridapDistributed
using PartitionedArrays
using MPI

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

# ## Discrete model
# To obtain a refined 2D parametric model, we pass $\ell$ levels of refinement:
ℓ = 2
radius = 1.0
omodel = ParametricOctreeDistributedDiscreteModel(ranks, radius; num_initial_uniform_refinements=ℓ)
model = omodel.parametric_dmodel

# ## Triangulation
# Now we extract the triangulated and the panel ids associated to each cell:
Ω = Triangulation(model)
panel_ids = get_panel_ids(model)


# ## FE Spaces
# Trial and test spaces relevant to compatible finite element pairings are defined
# using Gridap's high level API (refer to [Cotter et al. 2012](https://doi.org/10.1016/j.jcp.2012.05.020)).
order = 1
R = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,order+1); conformity=:H1)
H = TransientTrialFESpace(R)

V = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,order); conformity=:HDiv)
U = TransientTrialFESpace(V)

Q = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,order); conformity=:L2)
P = TransientTrialFESpace(Q)

# The multifield FE spaces for the prognostic and diagonstic variables are:
X_prog = MultiFieldFESpace([U,P]) # u, p
Y_prog = MultiFieldFESpace([V,Q]) # u, p

X_diag = MultiFieldFESpace([H,U,P]) # q, F, Φ
Y_diag = MultiFieldFESpace([R,V,Q]) # q, F, Φ

# ## Initial conditions
# The initial condition for velocity, fluid depth and Coriolis are defined
# as per the Williamson 2 test case, in [Williamson et al. 1992](https://doi.org/10.1016/S0021-9991(05)80016-6).
#
# The dimensional parameters are:
aₑ = 6.37e6 # m
gₑ = 9.8 # m/s^{-2}
ωₑ = 7.29e-5 # s^{-1}
Hₑ = 2.94e4/gₑ # m
uₑ = 2*π*aₑ/(12*24*3600) # m s^{-1}

# Consider a length scale of $a_{e}$ and time scale of $1/\omega$ yields the non-dimensional parameters:
gravity = gₑ*(1/ωₑ)^2/aₑ
ω = ωₑ*(1/ωₑ)
H₀ = Hₑ/aₑ
u₀ = uₑ/aₑ*(1/ωₑ)

# The initial conditions are defined as a function of the forward map as follows:
function u0(forward_map)
  function _u₀(α)
  ζ = 0.0
  xyz = forward_map(α)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  u     = u₀*(cos(ϕ)*cos(ζ) + cos(θ)*sin(ϕ)*sin(ζ))
  v     = - u₀*sin(θ)*sin(ζ)
  function _spherical_to_cartesian_matrix(θϕr)
    θ,ϕ,r = θϕr
    TensorValue(-sin(θ)       , cos(θ)       ,      0,
                -sin(ϕ)*cos(θ),-sin(ϕ)*sin(θ), cos(ϕ),
                 cos(ϕ)*cos(θ), cos(ϕ)*sin(θ), sin(ϕ))
  end
  _spherical_to_cartesian_matrix(θϕr)⋅VectorValue(u,v,0)
  end
end

function h0(forward_map)
  function _h₀(α)
  ζ = 0.0
  xyz = forward_map(α)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  h  = -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ)
  H₀- (ω*u₀  + 0.5*u₀*u₀)*h*h/gravity
  end
end

function f0(forward_map)
  function _f₀(α)
    ζ = 0.0
    xyz = forward_map(α)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    2.0*ω*( -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ) )
  end
end

# Then converted into a panelwise cellfield, where we extract the contravariant components
# for the velocity:
u_cf = panelwise_cellfield(piola(u0),Ω,panel_ids)
h_cf = panelwise_cellfield(h0,Ω,panel_ids)
f_cf = panelwise_cellfield(f0,Ω,panel_ids)


# ## Weak form
# To define the weak form, we require the metric and measure, as well as the
# the matrix that represents the perp operator.
# We use an increased degree of quadrature to exactly approximate the geometrical map included in the weak form.
g = panelwise_cellfield(metric,Ω,panel_ids)
ginv = panelwise_cellfield(inv_metric,Ω,panel_ids)
meas = panelwise_cellfield(sqrtg,Ω,panel_ids)
covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω,panel_ids)
Aperp = [0 -1
        1 0]
Rperp = TensorValue(Aperp)
Rperp_cf = CellField(Rperp,Ω)
dΩ = Measure(Ω,4*order)

# ### Diagnostic variables
# The weak forms for the diagnostic variables are:
resq(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( q*p*w*meas )dΩ - ∫( f_cf*w*meas )dΩ - ∫( (((Rperp_cf⋅u)⋅ginv)⋅∇(w))*meas)dΩ
resF(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( (F⋅(g⋅v))*(1/meas) )dΩ - ∫( p*(u⋅(g⋅v))*(1/meas))dΩ
resΦ(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( Φ*ψ*meas )dΩ - ∫( gravity*p*ψ*meas  )dΩ - ∫( 0.5*(u⋅(g⋅u))ψ*(1/meas))dΩ

res_y(t,((u,p),(q,F,Φ)),(w,v,ψ)) = resq(((u,p),(q,F,Φ)),(w,v,ψ)) + resF(((u,p),(q,F,Φ)),(w,v,ψ)) + resΦ(((u,p),(q,F,Φ)),(w,v,ψ))
jac_y(t,((u,p),(q,F,Φ)),(dq,dF,dΦ),(w,v,ψ)) = ∫( dq*p*w*meas )dΩ + ∫( (dF⋅(g⋅v))*(1/meas) )dΩ + ∫( dΦ*ψ*meas  )dΩ

# ### Prognostic variables
# The weak forms for the prognostic variables are as follows, where we use
# SUPG stabilisation in the velocity equation
res_p(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = ∫( r*(∇⋅F) )dΩ
res_u(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = ( ∫( ( q*( (Rperp_cf⋅F)⋅v)))dΩ
                              + ∫( -τ*( (q-q0)/dt)*( (Rperp_cf⋅F)⋅v))dΩ
                              + ∫( -τ*(u⋅∇(q))*( (Rperp_cf⋅ F)⋅v)*(1/meas))dΩ
                              - ∫( Φ*(∇⋅v) )dΩ
                  )

mass(t,(dut,dpt),(v,r)) = ∫( (dut⋅(g⋅v))*(1/meas) )dΩ + ∫( (dpt*r)*meas )dΩ
res_x(t,((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = res_u(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) + res_p(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0))
jac_x(t,((u,p),(q,F,Φ)),(du,dp),(v,r),(q0,F0,Φ0)) = ∫( -τ*(du⋅∇(q))*( (Rperp_cf⋅F)⋅v)*(1/meas))dΩ
jac_xt(t,((u,p),(q,F,Φ)),(dut,dpt),(v,r),(q0,F0,Φ0)) =  ∫( (dut⋅(g⋅v))*(1/meas) )dΩ + ∫((dpt*r)*meas )dΩ



# ## Transient problem
# The transient problem requires a linear solver for both the prognostic and diagnostic problem.
# In this example, we using a Conjugate Gradient method for both.
ls_diag = CGSolver(JacobiLinearSolver();rtol=1-12,verbose=i_am_main(ranks),name="diagnostic_solver")
ls_ode = CGSolver(JacobiLinearSolver();rtol=1-12,verbose=i_am_main(ranks),name="ode_solver")

# The transient parameters are:
t0 = 0.0
tF = 2*π
nsteps = 100
dt = tF/nsteps
τ = dt/2

# The transient operator, DAE operator and transient solution, integrated using SSPRK3, is defined as follows:
opT = TransientSemilinearFEOperator(mass,res_x,(jac_x,jac_xt),X_prog,Y_prog)
opFE = FEOperator(res_y,jac_y,X_diag,Y_diag)
opDAE = DAEFEOperator(opT,opFE,ls_diag)
ode_solver = RungeKutta(ls_ode,ls_ode,dt,:EXRK_SSP_3_3)
xh0 = interpolate_everywhere([u_cf,h_cf],X_prog)
solT  = solve(ode_solver,opDAE,t0,tF,xh0)


# ## Post processing
# The diagnostic solution can be extracted from the cache of the ode solution.
# To do so, we utilise the iterable functionality of ``solT`` as follows:

i_am_main(ranks) && mkpath("shallow_water_sol")
PartitionedArrays.barrier(ranks)

it = iterate(solT)
while !isnothing(it)
  data, state = it
  t, xh = data
  odeopcache = state[2][5][2]
  yh = odeopcache.diagnostics

  uh,ph = xh
  qh,Fh,Φh = yh

  i_am_main(ranks) && println(t)

  writevtk(Ω,"shallow_water_sol/solT_$t.vtu",
      cellfields=["vel"=>covariant_basis_cf⋅(1/meas*uh),"p"=>ph,"vort"=>qh],
      append=false,geo_map=latlon_geo_map_func(Ω))

  global it = iterate(solT, state)
end
