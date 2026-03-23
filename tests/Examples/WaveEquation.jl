# # Linear wave equation on the cubed sphere manifold
#
# This example solves the linear wave equation, given by
#
# ```math
# \begin{align*}
# \widetilde{\boldsymbol{u}} + \nabla_{\gamma} \varphi &= \widetilde{\boldsymbol{f}}_1 \quad \text{in} \quad \gamma, \\
# \widetilde{\varphi} + \nabla_{\gamma}\cdot \widetilde{\boldsymbol{u}}  &= \widetilde{f}_2 \quad \text{in} \quad \gamma,
# \end{align*}
# ```
#
# where $\gamma$ is the cubed sphere manifold, $\widetilde{\varphi}, \widetilde{f}_2: \gamma \rightarrow \mathbb{R}$
# are scalar valued functions defined in the ambient space of the manifold,
# $\widetilde{\boldsymbol{u}}, \widetilde{\boldsymbol{f}}_1 \in T_p \gamma $ are vector valued functions defined in the
# tangent space of the manifold,
# $\nabla_{\gamma}$ is the surface gradient operator, and
# $\nabla_{\gamma}\cdot$ is the surface divergence operator.
#
# We use $H(\mathrm{div})$ vector and $L^2$ scalar finite-elements, to solve in the parametric space of the cubed sphere.
# The weak formulation in the parametric space is: find $\boldsymbol{u} \in H(\mathrm{div},\mathcal{V})$, $\varphi \in L^2(\mathcal{V})$ such that
# ```math
# \begin{align*}
# \int_{\mathcal{V}} \boldsymbol{u}\cdot(g \boldsymbol{v}) \sqrt{g}
# - \int_{\mathcal{V}} \varphi \nabla\cdot(\sqrt{g} \boldsymbol{v})
#    &= \int_{\mathcal{V}} \boldsymbol{f}_1\cdot(g \boldsymbol{v}) \sqrt{g}
# \qquad   \forall \boldsymbol{v} \in H(\mathrm{div},\mathcal{V}) \\
# \int_{\mathcal{V}} \varphi \psi \sqrt{g}
# + \int_{\mathcal{V}} \psi \nabla\cdot(\sqrt{g} \boldsymbol{u})
#    &= \int_{\mathcal{V}} f_2 \psi \sqrt{g}
# \qquad   \forall \psi \in L^2(\mathcal{V})

# \end{align*}
# ```
# where $\mathcal{V}$ is the parametric space $f_1\in \mathbb{R}^n$, $f_2: \mathcal{V} \rightarrow \mathbb{R}$,
# $g$ is the Riemannian metric associated to the geometrical map $\sigma: \mathcal{V} \rightarrow \gamma$,
# and $\sqrt{g} = (\det{g})^{1/2}$ is the measure.


# ## Set up
# First load all required pacakges. In this example, we will use a distributed model, and the
# iterative solver provided by GridapSolvers.jl. So we also initialise MPI
using GridapGeosciences
using Gridap
using GridapP4est
using PartitionedArrays
using MPI

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

# ## Discrete model

# To obtain a refined 3D parametric model, we pass $\ell$ levels of refinement to the vertical and horizontal:
ℓ = 2
octree3_model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
                        num_horizontal_uniform_refinements=ℓ,
                        num_vertical_uniform_refinements=ℓ);
model = octree3_model.parametric_dmodel

# We can visualise the triangulation in the ambient space of the 3D sphere
# by passing a cellwise array of geometrical maps to writevtk:
Ω = Triangulation(model)
writevtk(Ω,"sphere_model",append=false,geo_map=geo_map_func(Ω))

# The 3D cubed sphere model has tags associated to the bottom, top and intermediate
# boundary cells. We can visualise each componeont of the model by passing the appropriate tag
# to the BoundaryTriangulation constructor
Γ_top = BoundaryTriangulation(model,tags=["top_boundary"])
Γ_bottom = BoundaryTriangulation(model,tags=["bottom_boundary"])
Γ_intermediate = BoundaryTriangulation(model,tags=["intermediate_boundary"])

writevtk(Γ_bottom,"boundary_bottom",append=false,geo_map=geo_map_func(get_panel_ids(Γ_bottom)))
writevtk(Γ_top,"boundary_top",append=false,geo_map=geo_map_func(get_panel_ids(Γ_top)))
writevtk(Γ_intermediate,"boundary_intermediate",append=false,geo_map=geo_map_func(get_panel_ids(Γ_intermediate)))



# ## FE Spaces
# Now that we have a discrete model, we define trial and test spaces using Gridap's high level API.
# We enforce zero dirichlet boundary conditions in the velocity space.
order = 1
tags = ["bottom_boundary",  "top_boundary"]

Q = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,order); conformity=:L2)
P = TrialFESpace(Q)
V = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,order); conformity=:HDiv, dirichlet_tags=tags)
U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

# The assoicated multifields are:
Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

# ## Manufactured solution
# We consider the method of manufactured solutions for analytic solutions defined in ambient space:
# ```math
# \begin{align*}
# \widetilde{\boldsymbol{u}} &= (-y,x,0)\\
# \widetilde{\varphi} &= \exp( -(y^2 + z^2) )
# \end{align*}
# ```
# These analytic solutions are defined as a function of the panel index $p$, as follows:

function u(p)
  function _u(α)
    x = ForwardMap(p)(α)
    VectorValue(-x[2],x[1],0.0)
  end
end

function φ(p)
  function _φ(α)
    x = ForwardMap(p)(α)
    exp(-(x[2]^2+x[3]^2))
  end
end

# Each cell is assigned a panel identifier, $p$, which is extracted as a cellwise array.
# This is used to generate a panelwise cellfield of the analytic functions:
panel_ids = get_panel_ids(model)
u_proj_cf = panelwise_cellfield(u,Ω,panel_ids)
phi_cf = panelwise_cellfield(φ,Ω,panel_ids)

# The cooresponding rhs forcing function is defined panelwise, as follows, where
# we extra the contravariant component for the velocity:
sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(u)),Ω,panel_ids)
sgrad_cf = panelwise_cellfield(sgrad(φ),Ω,panel_ids)
pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω,panel_ids)

f1 = pinvJ_cf ⋅ (u_proj_cf + sgrad_cf) # exact contravariant component
f2 = phi_cf + sdiv_cf


# ## Weak form
# The weak form is written as a mulitifield problem using  using Gridap's high level API.
# We use an increased degree of quadrature to exactly approximate the geometrical map included in the weak form.
meas = panelwise_cellfield(sqrtg,Ω,panel_ids)
_grad_meas = panelwise_cellfield(grad_meas,Ω,panel_ids)
g = panelwise_cellfield(metric,Ω,panel_ids)
dΩ = Measure(Ω,6*order)

biform_u((u,p),(v,q)) = ∫( (u⋅ (g⋅v))*meas )dΩ - ∫( p*(v⋅_grad_meas + meas*(∇⋅v) ) )dΩ
biform_p((u,p),(v,q)) = ∫( (p*q)*meas )dΩ + ∫( q*(u⋅_grad_meas + meas*(∇⋅u) )  )dΩ

biform((u,p),(v,q)) = biform_u((u,p),(v,q)) + biform_p((u,p),(v,q))
liform((v,q)) = ∫( f1⋅(g⋅v)*meas )dΩ + ∫( (f2*q)*meas )dΩ


# ## FE problem
# Now we can build the FE operator that represents the linearised wave equation.
op = AffineFEOperator(biform,liform,X,Y)

# We solve using a LU solver. One can also solve using MUMPS or another iterative solver.

ls = LUSolver()
A = get_matrix(op)
b = get_vector(op)
ns = numerical_setup(symbolic_setup(ls,A),A)
x = Gridap.Algebra.allocate_in_domain(A); fill!(x,0.0)
solve!(x,ns,b)
xh = FEFunction(X,x)
uh,ph = xh

# For the depth, the $L^2$ norm of the error between the exact and numerical soltuions is computed as
ep = phi_cf - ph
ep_l2 = sqrt(sum(∫((ep⋅ep)*meas)dΩ))

# For the velocity, the parametric velocity field is pushed to the ambient space
# via the covariant basis (i.e. Jacobian). Then the $L^2$ norm of the error between
# the exact and numerical soltuions is computed as
covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω,panel_ids)
uh_proj = covarient_basis_cf ⋅ uh
eu = u_proj_cf - uh_proj
eu_l2 = sqrt(sum(∫((eu⋅(g⋅eu))*meas)dΩ))


# ## Post processing
# The solution can be visualised in the ambient space by passing a
# cell-wise array of geometrical maps to Gridap's writevtk function
writevtk(Ω,"wave_equation",
        cellfields=["p"=>phi_cf,"ph"=>ph,"ep"=>ep,
                "uamb"=>u_proj_cf,"uamb_h"=>uh_proj, "eu"=>eu],
        append=false,geo_map=geo_map_func(Ω))
