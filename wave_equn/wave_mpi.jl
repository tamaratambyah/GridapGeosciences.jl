using GridapGeosciences
using PartitionedArrays
using MPI
using Gridap
using GridapDistributed
using DrWatson
using FillArrays
using GridapSolvers


nprocs = (1,1)
ranks = with_mpi() do distribute
  distribute(LinearIndices((prod(nprocs),)))
end

u_exact(x) = VectorValue(x[1]*(1-x[1]), x[2]*(1-x[2]),0.0)
p_exact(x) = 1.0 + 0.001*x[1]*(1-x[1]) + 0.001*x[2]*(1-x[2])
f_u(x) =  u_exact(x) + ∇(p_exact)(x)
f_p(x) = p_exact(x) + (∇⋅u_exact)(x)

p = 1
degree = 6

model = CubedSphereDiscreteModel(ranks,1;adaptive=true)
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
l2(e,dΩ) = sum(∫(e⊙e)dΩ)

reffe_u = ReferenceFE(raviart_thomas,Float64,p)
Vs = TestFESpace(model,reffe_u,conformity=:Hdiv)
Us = TrialFESpace(Vs)

reffe_p = ReferenceFE(lagrangian,Float64,p)
Qs = TestFESpace(model,reffe_p;conformity=:L2)
Ps = TrialFESpace(Qs)

tests = MultiFieldFESpace([Us,Ps])
trials = MultiFieldFESpace([Vs,Qs])


biform((u,p),(v,q),dΩ) = ( ∫( u⋅v - (∇⋅v)*p )dΩ
                         + ∫( p*q + (∇⋅u)*q )dΩ )
liform((v,q),dΩ) = ∫( f_u⋅v  + f_p⋅q )dΩ

a(u,v) = biform(u,v,dΩ)
l(v) = liform(v,dΩ)
op = AffineFEOperator(a,l,trials,tests)
A = get_matrix(op)
b = get_vector(op)


# smatrices, A, b = compute_hierarchy_matrices(trials,tests,biform,liform,degree)

solver = GMRESSolver(20;Pl=JacobiLinearSolver(),maxiter=5,atol=1e-14,rtol=1.e-14,verbose=true)
ns = numerical_setup(symbolic_setup(solver,A),A)

# Solve
x = pfill(0.0,partition(axes(A,2)))
solve!(x,ns,b)


# Error
Uh = get_fe_space(trials,1)
uh, ph = FEFunction(Uh,x)

l2(u_exact-uh,dΩ)
l2(p_exact-ph,dΩ)

writevtk(Ω,datadir("wave"),cellfields = ["uh"=>uh, "ph"=>ph,
                "u"=>CellField(u_exact,Ω), "p"=>CellField(p_exact,Ω)],append=false)
