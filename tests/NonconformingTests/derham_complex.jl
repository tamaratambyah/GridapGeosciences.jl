using GridapP4est
using Gridap
using MPI
using PartitionedArrays
using DrWatson

dir = datadir("derham")
!isdir(dir) && mkdir(dir)

us(x) = sin(x[1]*x[2])
uv(x) = VectorValue(sin(x[1]*x[2]),cos(x[1]*x[2]))

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))
coarse_model = CartesianDiscreteModel((0,1,0,1),(2,1))
model = OctreeDistributedDiscreteModel(ranks, coarse_model, 1)

ref_coarse_flags = map(_ -> [nothing_flag, nothing_flag, nothing_flag, nothing_flag, refine_flag, nothing_flag, refine_flag, nothing_flag], ranks)
model, _ = Gridap.Adaptivity.adapt(model, ref_coarse_flags);

p_fe = 0
lagreffe = ReferenceFE(lagrangian, Float64, p_fe+1)
rtreffe = ReferenceFE(raviart_thomas, Float64, p_fe)
dfreffe = ReferenceFE(lagrangian, Float64, p_fe)

V = FESpace(model,lagreffe; conformity=:H1)
U = TrialFESpace(V)

R = FESpace(model,rtreffe; conformity=:HDiv)
T = TrialFESpace(R)

D = FESpace(model,dfreffe; conformity=:L2)
E = TrialFESpace(D)

degree = 2*p_fe + 1
Ω = Triangulation(model)
dΩ = Measure(Ω,40)


# perp
function perp(vec)
  VectorValue(-vec[2],vec[1])
end

π0u = interpolate(us,U)
grap_perp_π0u = perp∘(∇(π0u))
π1u_grad_perp_u = interpolate(perp∘∇(us),R)
writevtk(Ω,dir*"/test_de_Rham1",cellfields=["grad_perp_π0u"=>grap_perp_π0u,
                                      "π1u_grad_perp_u"=>π1u_grad_perp_u,
                                      "error"=>π1u_grad_perp_u-grap_perp_π0u],append=false);

# a(u,v) = ∫( u⋅v )dΩ
# l(v) = ∫( uv⋅v )dΩ
# op = AffineFEOperator(a,l,T,R)
# π1u = solve(op)
π1u = interpolate(uv,T)
div_π1u = divergence(π1u)

e = π1u_grad_perp_u-grap_perp_π0u
sum(∫(e⋅e)dΩ)

div_u = divergence(uv)
a_l2(p,q) = ∫(p*q)dΩ
l_l2(q) = ∫(div_u*q)dΩ
op = AffineFEOperator(a_l2,l_l2,D,E)
π2_div_u = solve(op)
# π2_div_u = interpolate(div_u,E)
writevtk(Ω,dir*"/test_de_Rham2",nsubcells=10,cellfields=["div_π1u"=>div_π1u,
                                      "π2_div_u"=>π2_div_u,
                                      "error"=>div_π1u-π2_div_u],append=false);
e = div_π1u-π2_div_u
sum(∫( e⋅e )dΩ)
