using DrWatson
using Gridap
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using GridapGeosciences
using Test

using MPI
using PartitionedArrays
using GridapGeosciences.Distributed
using GridapP4est
using GridapDistributed


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

include("../Laplace/analytic_funcs.jl")
include("../convergence_tools.jl")
include("../missing_overloads.jl")

function _single_panel_model(model::DiscreteModelPortion)
  _grid = get_grid(model)
  ## construct the new grid by hand, so that specific cmaps are given
  nodes = collect1d(Geometry.get_node_coordinates(_grid)) # these are just junk nodes, never used
  ctype = collect1d(get_cell_type(_grid))
  cmaps = collect1d(get_cell_map(_grid))
  grid = Gridap.Geometry.UnstructuredGrid(nodes,get_cell_node_ids(_grid),
            get_reffes(_grid),ctype,OrientationStyle(_grid),
                      nothing,cmaps)

  topo = UnstructuredGridTopology(get_grid_topology(model))
  labels = FaceLabeling(topo)

  return UnstructuredDiscreteModel(grid,topo,labels)
end

function single_panel_laplace_beltrami_solver(ranks,
  nc,
  p_fe::Int,dir::String,f::Function,pid::Int,degree::Int,ls=LUSolver(),return_vtk=false)

  # ranks = get_ranks(panel_model)

  # nc = Int(sqrt(num_cells(panel_model)/6))
  # println(nc)

  panel_f = f(pid)
  panel_inv_metric = inv_metric(pid)
  panel_sqrtg = sqrtg(pid)
  panel_surflap = surflap(f)(pid)

  # 3D model
  # single_panel_model = CartesianDiscreteModel((0,1,-π/4,π/4,-π/4,π/4),(nc,nc,nc))

  # 2D model
  # panel_ids = get_panel_ids(panel_model)
  # mask = panel_ids .== pid
  # model = DiscreteModelPortion(panel_model,mask)
  # single_panel_model = _single_panel_model(model)
  single_panel_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(nc,nc))

  e = laplace_beltrami_manufactured(single_panel_model,p_fe,degree,dir,
        panel_f,panel_inv_metric,panel_sqrtg,panel_surflap,
        ls,return_vtk)

  ### convergence output for DrWatson
  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)
  n = nc
  output = @strdict e n p_fe degree pid
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("single_panel_nref$(n)_quad$(degree).jld2")), output)



  return e
end

function laplace_beltrami_manufactured(model,p_fe::Int,degree::Int,dir::String,
  f::Function,inv_metric::Function,sqrtg::Function,surflap::Function,
  ls=LUSolver(),return_vtk=false)

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1,
    dirichlet_tags="boundary")
  U = TrialFESpace(V,f)

  f_cf = CellField(f,Ω)
  inv_metric_cf = CellField(inv_metric,Ω)
  meas_cf = CellField(sqrtg,Ω)
  slap_panel_cf =  CellField(surflap,Ω)

  # i_am_main(ranks) && println("Zeromean: ", sum(∫(f_panel_cf*meas_cf)dΩ))
  sum(∫(f_cf*meas_cf)dΩ) < 1e-14

  rhs_cf = - slap_panel_cf

  poisson_biform(u,v) =  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
  poisson_liform(v) = ∫(  (rhs_cf*v)*meas_cf )dΩ
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

  # uh = solve(ls,op)

  ## for pvectors, the ghost may not be in the prange of the get_matrix
  ## This causes issues with GridapSolvers Krylov solvers, in the allocation of x
  ## To avoid, allocate x based on the domain of A
  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  uh = FEFunction(U,x)

  dΩ_error = Measure(Ω,16)
  e = l2(f_cf-uh,meas_cf,dΩ_error)

  if return_vtk
    cell_geo_map = geo_map_func(map(x->pid,collect(1:num_cells(model))))
    panel_cfs = [f_cf,uh,f_cf-uh,rhs_cf]
    labels = ["u","uh","eu","rhs"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω,dir*"/Single_panel",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return e
end

ranks=[true]
ns = [4]#collect(1:4)
ls = LUSolver()
p_fe = 2
pid = 1 # single panel id to consider

f(p::Int) = x-> (π/4 + x[1])*(π/4 - x[1]) + (π/4 + x[2])*(π/4 - x[2])

degrees = [4]#[i*p_fe for i=1:5]

tol = 1e-16
return_vtk = true

######### plot against n
dir = datadir("single_panel_2D")
!isdir(dir) && mkdir(dir)

plot()
for (i,degree) in enumerate(degrees)
  errors = Vector{Float64}(undef,length(ns))
  for (j,n) in enumerate(ns)
    errors[j] = single_panel_laplace_beltrami_solver(ranks,2^n,p_fe,dir,f,pid,degree,ls,return_vtk)
  end
  plot!(ns,errors,
    lw=3,markersize=6, markershape=:circle,
    label="quad degree = $degree")
  plot!(yscale=:log10,framestyle=:box,xlabel="ref lvl",ylabel="L2(u - uh)"
          )
end
plot!(ns,tol*ones(length(ns)),
      lw=2,ls=:dash,color=:black,label="tol = $tol")
plot!(show=true)
savefig(dir*"/quad_degree_3D")
