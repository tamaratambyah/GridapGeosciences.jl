"""
Solve the following on periodic meshes
u   + ∇p  = f
Cp  + ∇⋅u = g
where C = 1/c
"""

using DataStructures
using DrWatson
using DataFrames

using Test
using LinearAlgebra
using FillArrays, BlockArrays

using Gridap
using Gridap.ReferenceFEs, Gridap.Algebra, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.MultiField, Gridap.Algebra
using PartitionedArrays
using GridapDistributed

using GridapSolvers
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools, GridapSolvers.PatchBasedSmoothers
using GridapSolvers.BlockSolvers: LinearSystemBlock, BiformBlock, BlockTriangularSolver
using MPI
using GridapP4est

function adapt_model(ranks,model::OctreeDistributedDiscreteModel)
  cell_partition=get_cell_gids(model.dmodel)
  ref_flags=map(ranks,partition(cell_partition)) do rank,indices
      flags=zeros(Cint,length(indices))
      flags.=refine_flag
  end
  ref_model, adaptivity_glue = Gridap.Adaptivity.adapt(model,ref_flags)
  return ref_model, adaptivity_glue
end


function Gridap.CellData.get_triangulation(a::GridapDistributed.DistributedMultiFieldCellField)
  trians = map(get_triangulation,a.field_fe_fun)
  # @check all(map(t -> t === first(trians), trians))
  return first(trians)
end

# u_exact(x) = VectorValue(x[2]*(1-x[2]),x[1]*(1-x[1]))
# p_exact(x) = 1 + x[1]*(1-x[1])

u_exact(x) = VectorValue(sin(2*π*x[1])*cos(2*π*x[2]), -sin(2*π*x[2])*cos(2*π*x[1]) )
p_exact(x) = 1 + 0.1*sin(2*π*x[1])


function convergence(ranks;c,α,n,order,iters,itu,itp,dir,return_vtk,simName)

  dir_convergence = dir*"/convergence"
  !isdir(dir_convergence) && mkdir(dir_convergence)

  C = c ≠ 0 ? 1/c : 0.0

  coarse_model = CartesianDiscreteModel((0,1,0,1),(n,n),isperiodic=(true,true))
  dmodel = OctreeDistributedDiscreteModel(ranks, coarse_model)
  fmodel,glue = adapt_model(ranks,dmodel)

  mh = ModelHierarchy([fmodel.dmodel,dmodel.dmodel])
  model = get_model(mh,1)
  Ω = Triangulation(model)
  qdegree = 2*(order+2)
  dΩ = Measure(Ω,qdegree)

  V = TestFESpace(model,ReferenceFE(raviart_thomas,Float64,order);conformity=:Hdiv);
  U = TrialFESpace(V)
  Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:L2)

  X = MultiFieldFESpace([U,Q])
  Y = MultiFieldFESpace([V,Q])

  f(x) = u_exact(x) + ∇(p_exact)(x)
  g(x) = C*p_exact(x) + (∇⋅u_exact)(x)
  biform_u(u,v,dΩ) = ∫(u⋅v)dΩ
  biform_p(p,q,dΩ) = ∫(C*p*q )dΩ
  biform((u,p),(v,q),dΩ) = ( biform_u(u,v,dΩ) + biform_p(p,q,dΩ)
                            - ∫(divergence(v)*p)dΩ
                            + ∫(divergence(u)*q)dΩ
                            )
  liform((v,q),dΩ) = ∫(v⋅f)dΩ + ∫(q*g)dΩ

  a(u,v) = biform(u,v,dΩ)
  l(v) = liform(v,dΩ)
  op = AffineFEOperator(a,l,X,Y)
  A, b = get_matrix(op), get_vector(op);

  ##### Preconditioned external solver
  ls = GMRESSolver(40;Pr=JacobiLinearSolver(),Pl=nothing,maxiter=1000,rtol=1.e-8,verbose=i_am_main(ranks))
  ns = numerical_setup(symbolic_setup(ls,A),A)

  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)

  gmg_iters = Bool(itu) ? solver_u.log.num_iters : 0
  cg_iters = Bool(itp) ? solver_p.log.num_iters : 0
  kylov_iters = ls.log.num_iters

  xh = FEFunction(X,x)
  uh,ph = xh

  dΩ_error = Measure(Ω,2*qdegree)
  eu = uh-u_exact
  Eu = sqrt(sum(∫( eu⋅eu )dΩ_error )) # 1e-8

  ep = ph-p_exact
  Ep = sqrt(sum(∫( ep⋅ep )dΩ_error )) # O(1)

  if Bool(return_vtk)
    writevtk(Ω,dir*"/$(simName)",
      cellfields= ["uh"=>uh, "ph"=>ph, "eu"=>eu, "ep"=>ep, "u"=>u_exact],
      append=false)
  end

  output = @strdict Eu Ep gmg_iters cg_iters kylov_iters n order c α itu itp
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("$(simName).jld2")), output)

end
