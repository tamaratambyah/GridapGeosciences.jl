using Gridap
using GridapSolvers
using GridapGeosciences
using LinearAlgebra
using DrWatson

# p_exact(x) = 1 + x[1]x[2]x[3]
# u_exact(x) = VectorValue(
# x[2]x[3] - 3x[1]x[1]x[2]x[3] / (x[1]x[1] + x[2]x[2] + x[3]x[3]),
# x[1]x[3] - 3x[1]x[2]x[2]x[3] / (x[1]x[1] + x[2]x[2] + x[3]x[3]),
# x[1]x[2] - 3x[1]x[2]x[3]x[3] / (x[1]x[1] + x[2]x[2] + x[3]x[3])
# )


p_exact(x) = 1.0 +  0.1*x[2]x[2]
u_exact(x) = VectorValue(x[1]x[1] , x[2]x[2], 0.0 )

f_u(x) = u_exact(x) + ∇(p_exact)(x)
f_p(x) = p_exact(x) + (∇⋅u_exact)(x)


n = 4
p = 2
degree = 10#2*(p+1)

# model = CartesianDiscreteModel((0,1,0,1,0,1), (n,n,n))
model = CubedSphereDiscreteModel(n)
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)
dω = Measure(Ω,degree,ReferenceDomain())
l2(e,dΩ) = sum(∫(e⊙e)dΩ)

rt_reffe = ReferenceFE(raviart_thomas, Float64, p)
lg_reffe = ReferenceFE(lagrangian, Float64, p)

V = FESpace(model, rt_reffe; conformity=:Hdiv)
U = TrialFESpace(V)

Q = FESpace(model, lg_reffe; conformity=:L2)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U, P])
Y = MultiFieldFESpace([V, Q])

# project exact solution
a0((u,p),(v,q)) = ∫( u⋅v + p*q )dΩ
l0((v,q)) = ∫( u_exact⋅v + p_exact*q )dΩ
op = AffineFEOperator(a0,l0,X,Y)
u_ex, p_ex = solve(LUSolver(),op)

_f_u = u_ex + ∇(p_ex)
_f_p = p_ex + (∇⋅u_ex)

# evaluate(f_p_ex,get_cell_points(Ω))

# solve wave equation
a((u, p), (v, q)) = ∫( u⋅v)dΩ - ∫( (∇⋅v)* p)dΩ + ∫( p*q)dΩ + ∫( (∇⋅u)*q )dΩ
l((v, q)) = ∫( _f_u⋅v +  _f_p*q  )dΩ


_a((u, p), (v, q)) = ∫( u⋅v)dΩ - ∫( (DIV(v))* p)dω + ∫( p*q)dΩ + ∫( (DIV(u))*q )dω
_l((v, q)) = ∫( u_ex⋅v + ∇(p_ex)⋅v + p_ex*q  )dΩ + ∫( (DIV(u_ex))*q  )dω

op = AffineFEOperator(a, l, X, Y)

# PD = PatchDecomposition(model)
# P = GridapSolvers.VankaSolver(X,PD)

# P = GridapSolvers.VankaSolver(Y)
# solver = GMRESSolver(20;Pl=P,maxiter=1000,atol=1e-14,rtol=1.e-14,verbose=true)
solver = LUSolver()


A = get_matrix(op)
ns = numerical_setup(symbolic_setup(solver,A),A)
x = zeros(axes(A,2))
b = get_vector(op)
solve!(x,ns,b)

xh = solve(solver,op)

uh, ph = xh
l2(uh - u_ex,dΩ)
l2(ph - p_ex,dΩ)
writevtk(Ω,joinpath(datadir("models"),"cubed_sphere_wave"),
                cellfields=["uh"=>uh,"ph"=>ph,
                            "u_ex"=>u_ex, "p_ex"=>p_ex,
                            "u"=>CellField(u_exact,Ω),"p"=>CellField(p_exact,Ω)],append=false)
