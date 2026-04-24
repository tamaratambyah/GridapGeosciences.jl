""" Linear advection equation (flux form) with manufactured solution
u + ∇ᵧ⋅(βu) = f
Solve with dG upwinding as per Brezzi 2004 paper
Replicate test in Section 5.4 of Rognes2013 paper
"""

module AdvectionDGUpwinding

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using Test
using SparseArrays
using DataFrames
using CSV

include("advection_funcs.jl")
include("../convergence_tools.jl")

function my_mean( Bu_n::SkeletonPair)
  plus  = ( Bu_n.plus)
  minus = ( Bu_n.minus)
  0.5*( plus - minus  )
end

function _my_mean(j::SkeletonPair,vel::CellField,u::CellField)
  0.5*( (j.plus⋅vel.plus)*u.plus + (j.minus⋅vel.minus)*u.minus )
end

function _my_other_mean(j::SkeletonPair,vel::CellField,u::CellField,meas)
  0.5*( (j.plus⋅vel.plus)*u.plus*meas.plus + (j.minus⋅vel.minus)*u.minus*meas.minus )
end

function my_jump(j::SkeletonPair,ginv::SkeletonPair,n::SkeletonPair,u::CellField)
  u.plus*(j.plus⋅(ginv.plus⋅n.plus) ) + u.minus*(j.minus⋅(ginv.minus⋅n.minus) )
end


################################################################################
#### Steady with manufactured solutions
################################################################################
function advection_dg_solver(
  panel_model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}},
  p_fe::Int,dir::String,u::Function,vX::Function,uvX::Function,ls=LUSolver(),return_vtk=false)

  das = FullyAssembledRows()
  # das = SubAssembledRows()

  ranks = get_ranks(panel_model)

  i_am_main(ranks) && println("Assembly strategy: $das")

  lvl = nref(nc(panel_model))
  i_am_main(ranks) && println("nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  degree = 2*(p_fe + 1)

  Ω_panel = Triangulation(das,panel_model)
  dΩ = Measure(Ω_panel,degree)

  Λ = SkeletonTriangulation(das,panel_model)
  dΛ = Measure(Λ,degree)
  n_Λ = get_normal_vector(Λ)

  _rhs(p) = αβ -> u(p)(αβ) + surfdiv(contra_v(uvX))(p)(αβ)

  v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)
  rhs_cf = panelwise_cellfield(_rhs,Ω_panel,panel_ids)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  # hard code RT space as order 1 -- for velocity
  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
  U = TrialFESpace(V)


  vel = v_contr_cf

  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  meas_cf_skel = panelwise_cellfield(sqrtg,Λ)

  a_Ω(u,v) = ∫( (u*v)*meas_cf )dΩ - ∫( (u*(∇(v)⋅vel) )*meas_cf )dΩ

  ### volume stabilisation term
  a_s1(u,v) = ∫( my_mean((vel*u)⋅n_Λ)*jump(v)*meas_cf_skel.plus   )dΛ
  # a_s1(u,v) = ∫( my_mean((vel*u)⋅n_Λ)*jump(v)*_meas_cf   )dΛ

  # jac_cf = panelwise_cellfield(forward_jacobian,Λ)
  # ginv_cf = panelwise_cellfield(inv_metric,Λ)
  # a_s1(u,v) = ∫( _my_mean(jac_cf,vel,u)⋅my_jump(jac_cf,ginv_cf,n_Λ,v)*meas_cf_skel.plus   )dΛ
  # a_s1(u,v) = ∫( _my_other_mean(jac_cf,vel,u,meas_cf_skel)⋅my_jump(jac_cf,ginv_cf,n_Λ,v)   )dΛ

  ### upwinding stabilisation term
  upwind = abs( (vel⋅ n_Λ).plus)
  a_s2(u,v) = ∫(  0.5*(upwind)*jump(u)*jump(v)*meas_cf_skel.plus   )dΛ
  # a_s2(u,v) = ∫(  0.5*(upwind)*jump(u)*jump(v)*_meas_cf   )dΛ

  # cell_geo_map = geo_map_func(panel_ids)
  # n = pushforward_normal(Λ,cell_geo_map)
  # n = pushforward_normal(Λ)
  # a_s2(u,v) = ∫(  (0.5*(upwind)*jump(u*n)⋅jump(v*n))*meas_cf_skel.plus   )dΛ

  assem = SparseMatrixAssembler(P,Q,das)
  biform_advection(p,q) =  a_Ω(p,q) + a_s1(p,q) + a_s2(p,q)
  liform_advection(q) = ∫( (rhs_cf*q)*meas_cf )dΩ

  op = AffineFEOperator(biform_advection,liform_advection,P,Q,assem)

  # uh = solve(ls,op)
  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  uh = FEFunction(P,x)

  function save_sparse_matrix(M)
    I, J, V = SparseArrays.findnz(M)
    println("I = ", I)
    df = DataFrame([:I => I, :J => J, :V => V])
    CSV.write("M_nref$lvl.csv", df)
  end

  function save_vector(v,name)
    df = DataFrame(col=v)
    CSV.write(name, df)
  end

  # Atrivial = PartitionedArrays.to_trivial_partition(A)
  # btrivial = PartitionedArrays.to_trivial_partition(b, partition(axes(Atrivial,1)))

  # if (i_am_main(ranks))
  #   Aloc=partition(Atrivial).item_ref[]
  #   println("Matrix size: ",size(Aloc))
  #   save_sparse_matrix(Aloc)
  # end

  # if (i_am_main(ranks))
  #   bloc=partition(btrivial).item_ref[]
  #   save_vector(bloc,"b$(length(ranks)).csv")
  # end

  # if (i_am_main(ranks))
  #   Aloc=partition(btrivial).item_ref[]
  #   println("Matrix size: ",size(Aloc))
  #   save_sparse_vector(Aloc)
  # end

  Ωo = Triangulation(panel_model)
  dΩo = Measure(Ωo,degree)
  eu = l2((uh-u_cf),meas_cf,dΩo)
  i_am_main(ranks) && println("Error = ",eu)

  if return_vtk
    _Ω_panel = Triangulation(panel_model)
    lvl = nref(nc(panel_model))
    cell_geo_map = geo_map_func(get_panel_ids(_Ω_panel))
    labels = ["uh","u","eu"]
    panel_cfs = [uh,u_cf,uh-u_cf]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe", cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

    ### convergence output for DrWatson
    dir_convergence = dir*"/convergence"
    (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

    n = nc(panel_model)
    dxx = dx(panel_model)
    output = @strdict eu n dxx p_fe lvl
    i_am_main(ranks) && safesave(datadir(dir_convergence, ("advection_dg_nref$(lvl)_p$p_fe.jld2")), output)

  return eu,false,false
end

################################################################################
#### Auto convergence test
################################################################################
function main(distribute,nprocs;octree=false)
  ranks = distribute(LinearIndices((nprocs,)))
  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Auto conference test: AdvectionDGUpwinding")

  n_ref_lvls = 5
  ps = [1,2,3]
  ls = LUSolver()

  vX = panel_to_cartesian(tangent_vec(vecX))
  u = panel_to_cartesian(u0)
  uvX = panel_to_cartesian(u0vecX)

  dir = foldername("AdvectionDGConvergence",octree,false)
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  models = get_models(ranks,nprocs,n_ref_lvls;threedims=false,octree=octree)

  i_am_main(ranks) && println("advection_dg_convergence_func")
  p_convergence_test(ranks,ps,models,advection_dg_solver,dir,u,vX,uvX,ls,true)

  i_am_main(ranks) && println("--DONE--")
end



################################################################################
#### Convergence test with plots
################################################################################
function advection_dg_convergence_test(ranks,nprocs,u::Function,vX::Function,uvX::Function,
  n_ref_lvls=4,ps=[1],ls=LUSolver(),return_vtk=false)

  # serial models
  models  = get_refined_models(n_ref_lvls)

  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    models,  = get_distributed_refined_models(ranks,nprocs,models)
  end

  simName = "advection_dg_convergence_func"
  i_am_main(ranks) && println(simName)

  dir = datadir("AdvectionDG")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  errors = Vector{Vector{Float64}}(undef,length(ps))
  ns = Vector{Vector{Float64}}(undef,length(ps))
  dxs = Vector{Vector{Float64}}(undef,length(ps))
  slopes = Vector{Float64}(undef,length(ps))

  for (i,p_fe) in enumerate(ps)
    println("p_fe = $p_fe")
    errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,advection_dg_solver,p_fe,dir,u,vX,uvX,ls,return_vtk)
  end

  i_am_main(ranks) && print_convergence_results(errors,ns,dxs,slopes,ps)

  output = @strdict errors ns dxs slopes ps
  i_am_main(ranks) && safesave(datadir(dir, ("$simName.jld2")), output)

  i_am_main(ranks) && plot_convergence_from_saved(dir,simName)


end


end #module
