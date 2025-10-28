using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using Gridap.Algebra

using DrWatson
dir = datadir("Distributed")
!isdir(dir) && mkdir(dir)

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

num_horizontal_uniform_refinements = 2
num_vertical_uniform_refinements = 1
model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
	                                       num_horizontal_uniform_refinements=num_horizontal_uniform_refinements,
                                           num_vertical_uniform_refinements=num_vertical_uniform_refinements);


include("../convergence_tools.jl")
panel_model = model.parametric_dmodel
p_fe = 1
ls = LUSolver()
return_vtk = true
fXYZ(XYZ) =  XYZ[1]*XYZ[2]*XYZ[3]
f = panel_to_cartesian(fXYZ)

γαβ = Point(0.0,π/4,π/4)
αβ = map(x->Point(x[2],x[3]),γαβ)

surflap(f)(1)(γαβ)




panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,2*p_fe+1)

V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1, constraint=:zeromean)
U = TrialFESpace(V)

f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)

cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/laplace_beltrami",
cellfields=["u"=>f_panel_cf,"slapl"=>slap_panel_cf],
append=false,geo_map=cell_geo_map)


# i_am_main(ranks) && println("Zeromean: ", sum(∫(f_panel_cf*meas_cf)dΩ))
sum(∫(f_panel_cf*meas_cf)dΩ) < 1e-14
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

# l2(e,dΩ) = sum(∫( e⋅e )dΩ)
e = l2(f_panel_cf-uh,dΩ)

panel_cfs = [f_panel_cf,uh,f_panel_cf-uh]
labels = ["u","uh","eu"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/laplace_beltrami_sol",cellfields=cellfields,append=false,geo_map=cell_geo_map)


### convergence output for DrWatson
dir_convergence = dir*"/convergence"
(i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

n = nc(panel_model)
dxx = dx(nc(panel_model))
output = @strdict e n dxx p_fe lvl
i_am_main(ranks) && safesave(datadir(dir_convergence, ("laplace_beltrami_nref$(lvl)_p$p_fe.jld2")), output)

return e, false,false
