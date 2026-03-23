```@meta
EditURL = "../../../tests/Applications/LaplaceBeltrami.jl"
```

# Laplace Beltrami equation on the cubed sphere manifold

This example solves the Laplace Beltrami equation, given by

```math
\begin{align*}
-\Delta_{\gamma} \widetilde{u} &= \widetilde{f} \quad \text{in} \quad \gamma,
\end{align*}
```

where $\gamma$ is the cubed sphere manifold, $\widetilde{u}, \widetilde{f}: \gamma \rightarrow \mathbb{R}$
are scalar valued functions, defined in the ambient space of the manifold, and
$\Delta_{\gamma}$ is the Laplace Beltrami operator.

We use $H^1$ scalar finite-elements, to solve in the parametric space of the cubed sphere.
The weak formulation in the parametric space is: find $u \in H^1(\mathcal{V})$ such that
```math
\begin{align*}
\int_{\mathcal{V}} \mathbf{grad}~u \cdot (g^{-1} \mathbf{grad}~ v )~\sqrt{g} &= \int_{\mathcal{V}} f v ~\sqrt{g} \qquad   \forall v \in H^1(\mathcal{V})
\end{align*}
```
where $f: \mathcal{V} \rightarrow \mathbb{R}$, $g$ is the Riemannian metric, and $\sqrt{g} = (\det{g})^{1/2}$ is the measure.

## Set up
First load all required pacakges. In this example, we will use a serial model, and the
 basic LU solver provided in Gridap. One can use a distributed model and solvers.

````julia 
using GridapGeosciences
using Gridap
using GridapSolvers
````

## Discrete model

To obtain a refined parametric model, we first define the coarse parametric model, and
then apply $\ell$ levels of refinement, as follows:

````julia 
ℓ = 3
panel_model = coarse_parametric_model()
for n in collect(1:ℓ)
    panel_model = Gridap.Adaptivity.refine(panel_model)
end
````

Each cell is assigned a panel identifier, $p$. We extract these as a cellwise array:

````julia 
panel_ids = get_panel_ids(panel_model)
````

Using the panel ids, we can visualise the triangulation in the ambient space of the sphere
or in latitiude-longitude by passing a cellwise array of geometrical maps to writevtk:

````julia 
Ω_panel = Triangulation(panel_model)
writevtk(Ω_panel,"sphere_model",append=false,geo_map=geo_map_func(Ω_panel))
writevtk(Ω_panel,"latlon_model",append=false,geo_map=latlon_geo_map_func(Ω_panel))
````

## FE Spaces
Now that we have a discrete model, we define trial and test spaces using Gridap's high level API.
To remove the kernel, we use the zeromean constraint in the definition of the FE space.

````julia 
order = 1
V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,order); conformity=:H1, constraint=:zeromean)
U = TrialFESpace(V)
````

## Manufactured solution
We consider the method of manufactured solutions for analytic solution, $\widetilde{u} = xyz$.
This is defined as a function of the panel index $p$, as follows:

````julia 
function u(p)
  function _u(α)
    x = ForwardMap(p)(α)
    x[1]*x[2]*x[3]
  end
end
````

The cooresponding CellField and rhs forcing function is defined panelwise, as follows:

````julia 
f_panel_cf = panelwise_cellfield(u,Ω_panel,panel_ids)
slap_panel_cf =  panelwise_cellfield(surflap(u),Ω_panel,panel_ids)
rhs_cf = - slap_panel_cf
````

## Weak form
To define the weak form, we first obtain panelwise cellfields of the inverse metric and meas,
and then write the bilinear and linear forms using Gridap's high level API.
We use an increased degree of quadrature to exactly approximate the geometrical map included in the weak form.

````julia 
inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
dΩ = Measure(Ω_panel,6*order)
poisson_biform(u,v) =  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
poisson_liform(v) = ∫(  (rhs_cf*v)*meas_cf )dΩ
````

## FE problem
Now we can build the FE operator that represents the Lapalce Beltrami equation,
and solve using LU factorisation

````julia 
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)
````

The L2 norm of the error between the exact and numerical soltuions is computed in the
parametric space

````julia 
e = f_panel_cf-uh
el2 = sqrt(sum(∫( (e⋅e)*meas_cf )dΩ))
````

## Post processing
The solution can be visualised in the ambient space by passing a
cell-wise array of geometrical maps to Gridap's writevtk function

````julia 
writevtk(Ω_panel,"laplace_beltrami",cellfields=["u"=>f_panel_cf,"uh"=>uh,"eu"=>e],append=false,geo_map=geo_map_func(Ω_panel))
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

