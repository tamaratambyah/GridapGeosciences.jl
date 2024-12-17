using Gridap
using GridapSolvers
using GridapGeosciences
using LinearAlgebra
using DrWatson


out_dir = datadir("Hdiv")

nc = 10
order = 0
degree = 6

model = CubedSphereDiscreteModel(nc)
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)
l2(e,dΩ) = sum(∫(e⊙e)dΩ)

rt_reffe = ReferenceFE(raviart_thomas, Float64, order)

V = FESpace(model, rt_reffe; conformity=:Hdiv)
U = TrialFESpace(V)

u_exact(x) = VectorValue(1.0, 2.0, 3.0)

n = get_normal_vector(Ω)
# n_norm = Operation(n -> n⋅n)(n)
unit_normal = Operation(n -> n/(n⋅n))(n) # unit_normal

function decompose_normal_tangential_components(u,n)
  # n == unit normal vector
  u_n = (u⋅n)*n
  u_t = u - u_n
  return u_n, u_t
end

# cell field
u_cf = CellField(u_exact,Ω)
u_cf_n, u_cf_t = decompose_normal_tangential_components(u_cf,unit_normal)

# interpolate
u_int = interpolate(u_exact,U)
u_int_n, u_int_t = decompose_normal_tangential_components(u_int,unit_normal)

# project onto Hdiv
a(u,v) = ∫( u⋅v )dΩ
function l(f,dΩ)
 (v) -> ∫( f⋅v  )dΩ
end

op = AffineFEOperator(a,l(u_exact,dΩ),U,V)
u_proj = solve(LUSolver(),op)
u_proj_n,u_proj_t = decompose_normal_tangential_components(u_proj,unit_normal)

# op_n = AffineFEOperator(a,l(u_cf_n,dΩ),U,V)
# _u_proj_n = solve(LUSolver(),op_n)
# u_proj_n,_ = decompose_normal_tangential_components(_u_proj_n,unit_normal)

# op_t = AffineFEOperator(a,l(u_cf_t,dΩ),U,V)
# _u_proj_t = solve(LUSolver(),op_t)
# _,u_proj_t = decompose_normal_tangential_components(_u_proj_t,unit_normal)


writevtk(Ω,joinpath(out_dir,"normal_tangent_decompose_sphere"),
                cellfields=["u_cf"=>u_cf,"u_int"=>u_int,"u_proj"=>u_proj,
                            "u_cf_n" => u_cf_n,"u_cf_t"=> u_cf_t,
                            "u_int_n"=>u_int_n,"u_int_t"=>u_int_t,
                            "u_proj_n" => u_proj_n,"u_proj_t"=> u_proj_t,
                            "err_proj_t"=> u_cf_t-u_proj_t,
                            "err_int_t"=> u_cf_t-u_int_t,
                           ],append=false)



l2(u_cf_t-u_proj_t,dΩ)
l2(u_cf_t-u_int_t,dΩ)
l2(u_proj_t-u_int_t,dΩ)
