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

function boundary_refinement(dmodel::GridapDistributed.DistributedDiscreteModel)
  ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
    flags=zeros(Cint,length(indices))
    flags.=nothing_flag

    flags[1]=refine_flag
    flags[own_length(indices)]=refine_flag

    # To create some unbalance
    if (rank%2==0 && own_length(indices)>1)
        flags[div(own_length(indices),2)]=refine_flag
    end
    flags
  end
end

### refine in the middle of the domain [0,1]^2
function middle_refinement(dmodel::GridapDistributed.DistributedDiscreteModel)

  ref_coarse_flags=map(partition(get_cell_gids(dmodel.dmodel)), local_views(dmodel) ) do indices,lmodel
    flags=zeros(Cint,length(indices))
    flags.=nothing_flag

    cmap = get_cell_map(get_grid(lmodel))
    ref_points = get_cell_ref_coordinates(lmodel)
    coords = lazy_map(evaluate,cmap,ref_points)

    for (i,xy) in enumerate(coords)
      x = map(x->x[1],xy)
      y = map(x->x[2],xy)
      if any(x .> 0.4 ) && any(x .< 0.6 ) && any(y .> 0.4) && any(y .< 0.6)
          flags[i] = refine_flag
      end
    end
    flags
  end
  return ref_coarse_flags
end

function test_transfer_ops_and_redistribute(ranks,
                                            dmodel::GridapDistributed.DistributedDiscreteModel{Dc},
                                            order,
                                            lvl::Int,periodic=false) where Dc

  ref_coarse_flags = middle_refinement(dmodel)
  # ref_coarse_flags = boundary_refinement(dmodel)
  fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);

  # Solve coarse
  xH,XH = if periodic
    solve_darcy_periodic(dmodel,order,lvl)
  else
    solve_darcy(dmodel,order,lvl)
  end
   uH,pH=xH
  UH,PH=XH

  # Solve fine
  xh,Xh = if periodic
    solve_darcy_periodic(fmodel,order,lvl+1)
  else
    solve_darcy(fmodel,order,lvl+1)
  end
  uh,ph=xh
  Uh,Ph=Xh

  Ωh = Triangulation(fmodel)
  degree = 2*(order+1)
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
  degree = 2*(order+1)
  dΩH = Measure(ΩH,degree)

  dΩhH = Measure(ΩH,Ωh,2*order)
  aHp(u,v) = ∫(v⋅u)*dΩH
  lHp(v)   = ∫(v⋅uh)*dΩhH
  oph     = AffineFEOperator(aHp,lHp,UH,UH)
  uhH     = solve(oph)
  e       = uH - uhH
  el2     = sqrt(sum( ∫( e⋅e )*dΩH ))

end


function solve_darcy(model::GridapDistributed.DistributedDiscreteModel{Dc},order,lvl::Int) where {Dc}

  dirichlet_tags=[5,6]
  neumanntags = [7,8]

  V = FESpace(model,
              ReferenceFE(raviart_thomas,Float64,order),
              conformity=:Hdiv,
              dirichlet_tags=dirichlet_tags)

  Q = FESpace(model,
              ReferenceFE(lagrangian,Float64,order);
              conformity=:L2)

  U = TrialFESpace(V,u_ex)
  P = TrialFESpace(Q)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  trian = Triangulation(model)
  degree = 2*(order+1)
  dΩ = Measure(trian,degree)

  btrian = Boundary(model,tags=neumanntags)
  degree = 2*(order+1)
  dΓ = Measure(btrian,degree)
  nb = get_normal_vector(btrian)

  a((u,p),(v,q)) = ( ∫( u⋅v - p*(∇⋅v)   )dΩ
                  + ∫( p*q + (∇⋅u)*q )dΩ )
  b((v,q)) = ∫( u_rhs⋅v + p_rhs*q  )dΩ - ∫((v⋅nb)*p_ex )dΓ


  op = AffineFEOperator(a,b,X,Y)
  xh = solve(op)

  uh, ph = xh
  eu = u_ex - uh
  ep = p_ex - ph
  writevtk(trian,dir*"/sol_lvl$lvl",cellfields=["u"=>u_ex,"p"=>p_ex,"uh"=>uh,"ph"=>ph, "eu"=>uh-u_ex, "ep"=>ph-p_ex],append=false)


  l2(v) = sqrt(sum(∫(v⋅v)*dΩ))
  h1(v) = sqrt(sum(∫(v*v + ∇(v)⋅∇(v))*dΩ))

  eu_l2 = l2(eu)
  ep_l2 = l2(ep)

  i_am_main(ranks) && println("[L2 ERROR Non-periodic mesh: lvl $lvl] uh = ", eu_l2)
  i_am_main(ranks) && println("[L2 ERROR Non-periodic mesh: lvl $lvl] ph = ", ep_l2)

  tol = 1.0e-6
  @test eu_l2 < tol
  @test ep_l2 < tol

  xh,X
end


function solve_darcy_periodic(model::GridapDistributed.DistributedDiscreteModel{Dc},order,lvl::Int) where {Dc}

  Ω = Triangulation(model)
  degree = 2*(order+1)
  dΩ = Measure(Ω, degree)

  V = TestFESpace(Ω,ReferenceFE(raviart_thomas,Float64,order),conformity=:HDiv)
  U = TrialFESpace(V)

  W = TestFESpace(Ω,ReferenceFE(lagrangian,Float64,order),conformity=:L2)
  R = TrialFESpace(W)

  X = MultiFieldFESpace([U,R])
  Y = MultiFieldFESpace([V,W])

  # Weak formulation for u
  # a((u,p),(v,q)) = ( ∫( u⋅v - p*(∇⋅v)   )dΩ
  #                   + ∫( p*q + (∇⋅u)*q )dΩ )
  # b((v,q)) = ∫( u_rhs⋅v + p_rhs*q  )dΩ
    a((u,p),(v,q)) = ( ∫( u⋅v   )dΩ
                    + ∫( p*q  )dΩ )
  b((v,q)) = ∫( u_ex⋅v + p_ex*q  )dΩ

  op = AffineFEOperator(a,b,X,Y)
  xh = solve(op)
  uh, ph = xh

  writevtk(Ω,dir*"/sol_periodic_lvl$lvl",cellfields=["u"=>u_ex,"p"=>p_ex,"uh"=>uh,"ph"=>ph, "eu"=>uh-u_ex, "ep"=>ph-p_ex],append=false)

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

# order = 1
# u_ex(x) = VectorValue(2*x[1],x[1]+x[2])
# p_ex(x) = x[1]-x[2]

order = 2
u_ex(x) = VectorValue(x[1]*(1-x[1]),0.0)
p_ex(x) = 1.0 + x[1]*(1-x[1])

u_rhs(x) = u_ex(x) + ∇(p_ex)(x)
p_rhs(x) = p_ex(x) + (∇⋅u_ex)(x)

### Non-periodic mesh
coarse_model = CartesianDiscreteModel((0,1,0,1),(3,3))
dmodel = OctreeDistributedDiscreteModel(ranks,coarse_model)#,2
test_transfer_ops_and_redistribute(ranks,dmodel,order,0,false)

### Periodic mesh
coarse_model = CartesianDiscreteModel((0,1,0,1),(5,5),isperiodic=(true,true))
dmodel = OctreeDistributedDiscreteModel(ranks,coarse_model)#,2)
test_transfer_ops_and_redistribute(ranks,dmodel,order,0,true)
