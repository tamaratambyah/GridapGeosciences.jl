using DrWatson
using Gridap

dir = datadir("Nonconfirming")
!isdir(dir) && mkdir(dir)

order = 2
# u_ex(x) = VectorValue(x[1]*(1-x[1]),x[2]*(1-x[2]))
# # p_ex(x) = 1.0 + x[1]*(1-x[1])
# function p_ex(x)
#   if x[1] < 0.5
#     return x[1]*(0.5-x[1])
#   else
#     return (x[1]-0.5)*(x[1]-1)
#   end
# end

# u_rhs(x) = u_ex(x) + ∇(p_ex)(x)
# p_rhs(x) = (∇⋅u_ex)(x)

model = CartesianDiscreteModel((0,1,0,1),(6,6),isperiodic=(true,true))

Ω = Triangulation(model)
degree = 2*(order+1)
dΩ = Measure(Ω, degree)

l2(u) = sqrt(sum(∫(u ⊙ u) * dΩ))

# sum(∫( p_rhs )dΩ)

V = TestFESpace(Ω,ReferenceFE(raviart_thomas,Float64,order),conformity=:HDiv)
U = TrialFESpace(V)

W0 = TestFESpace(Ω,ReferenceFE(lagrangian,Float64,order),conformity=:L2)
W = W0# Gridap.FESpaces.ZeroMeanFESpace(W0,dΩ)
R = TrialFESpace(W)

X = MultiFieldFESpace([U,R])
Y = MultiFieldFESpace([V,W])

# Weak formulation for u
biform((u,p),(v,q)) = ( ∫( u⋅v - p*(∇⋅v)   )dΩ
                  + ∫(  (∇⋅u)*q )dΩ )
# liform((v,q)) = ∫( u_rhs⋅v  )dΩ  + ∫( q*p_rhs )dΩ
# f(x) = VectorValue(x[1]*x[2],0.0)
f(x) = VectorValue(cos(x[1]*x[2]),0.0)
sum(∫(∇⋅f)dΩ)
liform((v,q)) = ∫( f⋅v  )dΩ


op = AffineFEOperator(biform,liform,X,Y)
A = get_matrix(op)
b = get_vector(op)
ns = numerical_setup(symbolic_setup(LUSolver(),A),A)
x = Gridap.Algebra.allocate_in_domain(A); fill!(x,0.0)
solve!(x,ns,b)
xh = FEFunction(X,x)
uh, ph = xh
sum(∫( ∇⋅uh )dΩ)
sum(∫( ph )dΩ)

# Error
# eu = u_ex - uh
# ep = p_ex - ph

# eu_l2 = l2(eu)
# ep_l2 = l2(ep)

writevtk(Ω,dir*"/mass_conservation",cellfields=["uh"=>uh,"ph"=>ph,"prhs"=>f, "divu"=>∇⋅uh],
append=false)


###
using GridapP4est
using PartitionedArrays
using GridapDistributed
using MPI
using Gridap.FESpaces

include("refinement_helpers.jl")

MPI.Init()
np = MPI.Comm_size(MPI.COMM_WORLD)
ranks = distribute_with_mpi(LinearIndices((np,)))

coarse_model = CartesianDiscreteModel((0,1,0,1),(6,6),isperiodic=(true,true))
dmodel = OctreeDistributedDiscreteModel(ranks,coarse_model)#,2)

ref_coarse_flags = middle_refinement(dmodel)
# ref_coarse_flags = boundary_refinement(dmodel)
fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);

model = fmodel
Ω = Triangulation(model)
degree = 2*(order+1)
dΩ = Measure(Ω, degree)

l2(u) = sqrt(sum(∫(u ⊙ u) * dΩ))

V = TestFESpace(Ω,ReferenceFE(raviart_thomas,Float64,order),conformity=:HDiv)
U = TrialFESpace(V)

W0 = TestFESpace(Ω,ReferenceFE(lagrangian,Float64,order),conformity=:L2)
W = W0# Gridap.FESpaces.ZeroMeanFESpace(W0,dΩ)
R = TrialFESpace(W)

X = MultiFieldFESpace([U,R])
Y = MultiFieldFESpace([V,W])

op = AffineFEOperator(biform,liform,X,Y)
A = get_matrix(op)
b = get_vector(op)
ns = numerical_setup(symbolic_setup(LUSolver(),A),A)
x = Gridap.Algebra.allocate_in_domain(A); fill!(x,0.0)
solve!(x,ns,b)

xh = FEFunction(X,x)
uh,ph = xh
sum(∫( ∇⋅uh )dΩ)
sum(∫( ph )dΩ)
sum(∫( ∇⋅f )dΩ )

writevtk(Ω,dir*"/mass_conservation_nonconforming",cellfields=["uh"=>uh,"ph"=>ph,"prhs"=>f,"divu"=>∇⋅uh],
append=false)
