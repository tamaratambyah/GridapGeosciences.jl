using GridapP4est
using Gridap
using PartitionedArrays
using GridapDistributed
using MPI
using Gridap.FESpaces
using FillArrays
using Test
using DrWatson

dir = datadir("Nonconfirming")
!isdir(dir) && mkdir(dir)
include("refinement_helpers.jl")

function test_transfer_ops_and_redistribute(ranks,
                                            dmodel::GridapDistributed.DistributedDiscreteModel{Dc},
                                            p_fe,
                                            lvl::Int) where Dc

  ref_coarse_flags = middle_refinement(dmodel)
  # ref_coarse_flags = boundary_refinement(dmodel)
  fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);

  # Solve coarse
  xH,XH = solve_darcy_periodic(dmodel,p_fe,lvl)

  uH,pH=xH
  UH,PH=XH

  # Solve fine
  xh,Xh = solve_darcy_periodic(fmodel,p_fe,lvl+1)
  uh,ph=xh
  Uh,Ph=Xh

  Ωh = Triangulation(fmodel)
  degree = 2*(p_fe+1)
  dΩh = Measure(Ωh,degree)

  # prolongation via interpolation
  uHh=interpolate(uH,Uh)
  e = uh - uHh
  el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
  tol=1e-6
  i_am_main(ranks) && println("[INTERPOLATION] el2 < tol: $(el2) < $(tol)")
  @assert el2 < tol

  # prolongation via L2-projection
  # Coarse FEFunction -> Fine FEFunction, by projection
  ahp(u,v)  = ∫(v⋅u)*dΩh
  lhp(v)    = ∫(v⋅uH)*dΩh
  oph      = AffineFEOperator(ahp,lhp,Uh,Uh)
  uHh      = solve(oph)
  e = uh - uHh
  el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
  i_am_main(ranks) && println("[L2 PROJECTION] el2 < tol: $(el2) < $(tol)")
  @assert el2 < tol

  # restriction via interpolation
  uhH=interpolate(uh,UH)
  e = uH - uhH
  el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
  i_am_main(ranks) && println("[INTERPOLATION] el2 < tol: $(el2) < $(tol)")
  @assert el2 < tol

  # restriction via L2-projection
  ΩH = Triangulation(dmodel)
  degree = 2*(p_fe+1)
  dΩH = Measure(ΩH,degree)

  dΩhH = Measure(ΩH,Ωh,2*p_fe)
  aHp(u,v) = ∫(v⋅u)*dΩH
  lHp(v)   = ∫(v⋅uh)*dΩhH
  oph     = AffineFEOperator(aHp,lHp,UH,UH)
  uhH     = solve(oph)
  e       = uH - uhH
  el2     = sqrt(sum( ∫( e⋅e )*dΩH ))

end


function solve_darcy_periodic(model::GridapDistributed.DistributedDiscreteModel{Dc},p_fe,lvl::Int) where {Dc}

  Ω = Triangulation(model)
  degree = 2*(p_fe+1)
  dΩ = Measure(Ω, degree)

  V = TestFESpace(Ω,ReferenceFE(raviart_thomas,Float64,p_fe),conformity=:HDiv)
  U = TrialFESpace(V)

  W = TestFESpace(Ω,ReferenceFE(lagrangian,Float64,p_fe),conformity=:L2)
  R = TrialFESpace(W)

  X = MultiFieldFESpace([U,R])
  Y = MultiFieldFESpace([V,W])

  # Weak formulation for u
  a((u,p),(v,q)) = ( ∫( u⋅v - p*(∇⋅v)   )dΩ
                    + ∫( p*q + (∇⋅u)*q )dΩ )
  b((v,q)) = ∫( u_rhs⋅v + p_rhs*q  )dΩ

  op = AffineFEOperator(a,b,X,Y)
  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(LUSolver(),A),A)
  x = Gridap.Algebra.allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)

  xh = solve(op)
  uh, ph = xh

  writevtk(Ω,dir*"/sol_darcy_periodic_lvl$lvl",
  cellfields=["u"=>u_ex,"p"=>p_ex,"uh"=>uh,"ph"=>ph, "eu"=>uh-u_ex, "ep"=>ph-p_ex],
  append=false,nsubcells=4)

  # Error
  eu = u_ex - uh
  ep = p_ex - ph
  l2(u) = sqrt(sum(∫(u ⊙ u) * dΩ))

  eu_l2 = l2(eu)
  ep_l2 = l2(ep)

  i_am_main(ranks) && println("[L2 ERROR Periodic mesh: lvl $lvl] uh = ", eu_l2)
  i_am_main(ranks) && println("[L2 ERROR Periodic mesh: lvl $lvl] ph = ", ep_l2)

  tol = 1.0e-6
  @test eu_l2 < tol
  @test ep_l2 < tol

  xh,X
end

MPI.Init()
np = MPI.Comm_size(MPI.COMM_WORLD)
ranks = distribute_with_mpi(LinearIndices((np,)))

p_fe = 2
u_ex(x) = VectorValue(x[1]*(1-x[1]),0.0)
p_ex(x) = 1.0 + x[1]*(1-x[1])

u_rhs(x) = u_ex(x) + ∇(p_ex)(x)
p_rhs(x) = p_ex(x) + (∇⋅u_ex)(x)

# ### Periodic mesh
coarse_model = CartesianDiscreteModel((0,1,0,1),(3,3),isperiodic=(true,true))
dmodel = OctreeDistributedDiscreteModel(ranks,coarse_model)#,2)

# ref_coarse_flags = middle_refinement(dmodel)
ref_coarse_flags = boundary_refinement(dmodel)
fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);

xH,XH = solve_darcy_periodic(dmodel,p_fe,0)
xh,Xh = solve_darcy_periodic(fmodel,p_fe,1)

test_transfer_ops_and_redistribute(ranks,dmodel,p_fe,0)
