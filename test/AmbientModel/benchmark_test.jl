"""
In this test, benchmark the assembly of the coriolis term in the parametric model
vs. the ambient model

Output: bm_ambient: 5.5 ms, mem 7.71 MiB, alloc 91262
        bm_panel: 1.63 ms,  mem 3.5 MiB,  alloc 25143
Conclusion: Parametric model is much better for performance for coriolis
most likely due to cross product.
"""

using GridapGeosciences
using Gridap
using BenchmarkTools



p_fe = 1
radius = 1



function get_setup(model,p_fe)
  degree = 5*(p_fe+1)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  Q = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  return Ω,dΩ,X,Y
end

function fcor_X(xyz)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    sin(ϕ)
  end

################################################################################
########## Ambient model
################################################################################
ambient_model = CubedSphereAmbientDiscreteModel(radius; num_initial_uniform_refinements=1)

Ω_ambient,dΩ_ambient,X_ambient,Y_ambient = get_setup(ambient_model,p_fe)
cor_cf_ambient = CellField(fcor_X,Ω_ambient)

## Here we construct the coriolis term on the surface: ∫( ̃f ( ̃k × ̃u  )  )dΩ
n_surf = get_surface_normal(Ω_ambient)
coriolis_term_ambient((u,p),(v,q)) = ∫( cor_cf_ambient*( ( n_surf × u)⋅v)  )dΩ_ambient

bm_ambient() = assemble_matrix(coriolis_term_ambient,X_ambient,Y_ambient)



################################################################################
########## Parametric model
################################################################################
panel_model = get_parametric_model(ambient_model)

Ω_panel,dΩ_panel,X_panel,Y_panel = get_setup(panel_model,p_fe)

# metric information
metric_cf = ParametricCellField(metric,Ω_panel)
meas_cf = ParametricCellField(sqrtg,Ω_panel)
cor_cf_panel = ParametricCellField(panel_to_cartesian(fcor_X),Ω_panel)



## Here we construct the coriolis term in 2D using the rotation matrx
function vecPerp(u)
  # u   = (u1, u2)
  # u^T = (-u2, u1)
  VectorValue(-u[2],u[1])
end
Aperp = [0 -1
        1 0]
Rperp = TensorValue(Aperp)
Rperp_cf = CellField(Rperp,Ω_panel)

coriolis_term_panel((u,p),(v,q)) = ∫( ( cor_cf_panel*( (Rperp_cf⋅ u)⋅v))  )dΩ_panel

bm_panel() = assemble_matrix(coriolis_term_panel,X_panel,Y_panel)


@benchmark bm_ambient()
@benchmark bm_panel()
