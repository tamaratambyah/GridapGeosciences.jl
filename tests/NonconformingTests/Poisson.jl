using MPI
using PartitionedArrays
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson
using Test

dir = datadir("Nonconfirming")
!isdir(dir) && mkdir(dir)

MPI.Init()
np = MPI.Comm_size(MPI.COMM_WORLD)
ranks = distribute_with_mpi(LinearIndices((np,)))

order = 2
function u_ex(x)
  if x[1] < 0.5
    return x[1]*(0.5-x[1])
  else
    return (x[1]-0.5)*(x[1]-1)
  end
end
f_ex(x) = u_ex(x) + Δ(u_ex)(x)
sum(∫( f_ex)dΩ  ) # check compatibility

coarse_model = CartesianDiscreteModel((0,1,0,1),(4,4),isperiodic=(true,true))
dmodel = OctreeDistributedDiscreteModel(ranks,coarse_model,2)

ref_coarse_flags=map(partition(get_cell_gids(dmodel.dmodel)), local_views(dmodel) ) do indices,lmodel
    flags=zeros(Cint,length(indices))
    flags.=nothing_flag

    cmap = get_cell_map(get_grid(lmodel))
    ref_points = get_cell_ref_coordinates(lmodel)
    coords = lazy_map(evaluate,cmap,ref_points)

    for (i,xy) in enumerate(coords)
      x = map(x->x[1],xy)
      y = map(x->x[2],xy)
      if any(x .> 0.3 ) && any(x .< 0.7 ) && any(y .> 0.3) && any(y .< 0.7)
          flags[i] = refine_flag
      end
    end
    flags
end
fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);
writevtk(Triangulation(fmodel),dir*"/fmodel",append=false)


model = fmodel.dmodel
# model = coarse_model

Ω = Triangulation(model)
degree = 2*(order+1)
dΩ = Measure(Ω,degree)


V = FESpace(Ω,ReferenceFE(lagrangian,Float64,order);conformity=:H1)
U = TrialFESpace(V)

a(u,v) = ∫( u*v )dΩ - ∫( ∇(u)⋅∇(v)  )dΩ
b(v) = ∫( v*f_ex )dΩ

op = AffineFEOperator(a,b,U,V)
uh = solve(op)


l2(v) = sqrt(sum(∫(v⋅v)*dΩ))
eu = u_ex - uh
eu_l2 = l2(eu)

writevtk(Ω,dir*"/sol",cellfields=["u"=>u_ex,"uh"=>uh,"eu"=>uh-u_ex],append=false)


########### dual form
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,order); conformity=:L2)
U = TrialFESpace(V)

T = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,order); conformity=:Hdiv)
S = TrialFESpace(T)

X = MultiFieldFESpace([S,U])
Y = MultiFieldFESpace([T,V])


biformX((s,u),(t,v)) = (  ∫( s⋅t)dΩ + ∫( (∇⋅t)*u )dΩ
                        + ∫( u*v )dΩ   + ∫( (∇⋅s)*v  )dΩ
                      )
liformY((t,v)) = ∫( f_ex*v )dΩ

op = AffineFEOperator(biformX,liformY,X,Y)
sh,uh = solve(LUSolver(),op)

e = uh - u_ex
l2(e)
writevtk(Ω,dir*"/sol",cellfields=["u"=>u_ex,"uh"=>uh,"eu"=>uh-u_ex],append=false)
