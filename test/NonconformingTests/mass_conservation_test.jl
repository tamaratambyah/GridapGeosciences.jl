using GridapP4est
# add Gridap#master
using Gridap
using MPI
using PartitionedArrays
using DrWatson
using GridapDistributed

include("refinement_helpers.jl")
include("../convergence_tools.jl")


dir = datadir("derham_trig")
!isdir(dir) && mkdir(dir)


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))
coarse_model = CartesianDiscreteModel((0,1,0,1),(4,4),isperiodic=(true,true))
dmodel = OctreeDistributedDiscreteModel(ranks, coarse_model)
ref_coarse_flags = initial_unbalance(dmodel)
fmodel, =Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);


# model = dmodel.dmodel.models.item
model = fmodel.dmodel.models.item

p_fe = 1
Ω = Triangulation(model)
dΩ = Measure(Ω,4)

Λ = SkeletonTriangulation(Ω)
dΛ = Measure(Λ,4)
n_Λ = get_normal_vector(Λ)


p_scalar(x) = 1 + 0.1*x[1]*(1-x[1])
u_vector(x) = VectorValue(x[1]*(1-x[1]), x[2]*(1-x[2]))
# u_vector(x) = VectorValue(sin(2*π*x[1])*cos(2*π*x[2]), -sin(2*π*x[2])*cos(2*π*x[1]))
div_u(x) = (∇⋅u_vector)(x)

V = FESpace(model, ReferenceFE(raviart_thomas, Float64, p_fe); conformity=:HDiv)
U = TrialFESpace(V)

W = FESpace(model, ReferenceFE(lagrangian, Float64, p_fe); conformity=:L2)
R = TrialFESpace(W)


u_h = interpolate(u_vector,U)

u_skel = CellField(u_vector,Λ)

jj = jump(u_skel)

writevtk(Ω,dir*"/velocity",
cellfields=["u"=>u_vector, "div_u"=>div_u, "uh"=>u_h,"p"=>p_scalar],append=false);

dd = u_h⋅n_Λ
ff = 0.5*(dd.plus+dd.minus)
writevtk(Λ,dir*"/velocity_jump",
      cellfields=["avF"=>mean(u_h),
      "plus"=>dd.plus,"minus"=>dd.minus, "avv"=>ff,
                  "jumpFn"=>jump(u_h⋅n_Λ),
                  "Fplus"=>((u_h⋅n_Λ).plus),
                  "Fminus"=>( (u_h⋅n_Λ).minus ),
                  "Fplusminus"=>((u_h⋅n_Λ).plus+(u_h⋅n_Λ).minus ),
                  "jumpplus"=>jump( (u_h⋅n_Λ).plus),
                  "jumpminus"=>jump( (u_h⋅n_Λ).minus) ],append=false)


### Hdiv conforming
a(p,q) = ∫( p*q )dΩ
a_div(q) = ∫( -1.0*divergence(u_h)*q )dΩ
b(q) =  a_div(q)
op = AffineFEOperator(a,b,R,W)
ph_hdiv = solve(LUSolver(),op)
sum(∫( ph )dΩ)



### in L2 with jump terms
rhs(x) = p_scalar(x) + (∇⋅u_vector)(x)
a(p,q) = ∫( p*q )dΩ - ∫( (gradient(q)⋅u_h)*p )dΩ + ∫( mean(u_h*p)⋅jump(q*n_Λ)  )dΛ
b(q) =  ∫( rhs*q )dΩ
op = AffineFEOperator(a,b,R,W)
ph = solve(LUSolver(),op)
sum(∫( ph )dΩ)

writevtk(Ω,dir*"/mass",
cellfields=["p"=>p_scalar, "ph_div"=>ph_hdiv, "ph"=>ph, "e"=>ph_hdiv-ph],append=false);


sum(a_skel(ph))
sum(∫( divergence(u_h)*ph )dΩ)
