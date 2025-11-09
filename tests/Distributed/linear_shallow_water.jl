using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using Gridap.Algebra

using DrWatson
dir = datadir("DistributedLinearisedSW")
!isdir(dir) && mkdir(dir)

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

include("../convergence_tools.jl")
include("williamson_funcs_3D.jl")

models = get_3D_octree_refined_models(ranks,3)
panel_model = models[1]
# panel_model = models[1]
p_fe = 1
ls = LUSolver()
return_vtk = true
ζ = 0.0

h = panel_to_cartesian(h₀(ζ))
vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
f = panel_to_cartesian(f₀(ζ))


# p_convergence_test(ranks,[p_fe],models,wave_solver,dir,h,vX,ls,true)
errors,ns,dxs,slopes = h_convergence_test(models,linear_shallow_water,p_fe,dir,h,vX,f,ls,true)

plot()
plot_convergence(errors,ns,dxs,slopes;leginf=["u","p"],colors=[palette(:tab10)[p_fe],palette(:tab10)[p_fe] ] )
savefig(dir*"/convergence_linear_shallow_wave_3D")



function linear_shallow_water(panel_model,p_fe::Int,dir::String,h::Function,vX::Function,f::Function,ls=LUSolver(),return_vtk=false)
Dc = num_cell_dims(panel_model)
nref(panel_model)
println("Dc = $Dc")

panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,2*(p_fe+1))

tags = ["bottom_boundary",  "top_boundary"]

Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)


V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])




covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω_panel,panel_ids)

h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)
cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)



# extract compoents 1 or 2 of contravariat vector, and construct contravariat components of vec perp
contra_v_comp3D(vecX::Function,p::Int,comp::Int) = αβ -> (forward_pinv_jacobian(p)(αβ)⋅ vecX(p)(αβ))[comp]
contra_v_comp3D(vecX::Function,comp::Int) = p -> contra_v_comp3D(vecX,p,comp)

contra_v_perp3D(vecX::Function,p::Int) = αβ -> sqrtg(p,αβ)*(
        inv_metric(p,αβ) ⋅ VectorValue(0.0, -contra_v_comp3D(vecX,p,3)(αβ), contra_v_comp3D(vecX,p,2)(αβ) ) )
contra_v_perp3D(vecX::Function) = p -> contra_v_perp3D(vecX,p)

u_perp_contra = panelwise_cellfield(contra_v_perp3D(vX),Ω_panel,panel_ids)
u_perp = covarient_basis_cf ⋅ u_perp_contra



sgrad_cf = panelwise_cellfield(sgrad(h),Ω_panel,panel_ids)
sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(vX)),Ω_panel,panel_ids)


# manufacture rhs functions
rhs_scalar = h_cf + sdiv_cf
rhs_vector = u_proj_cf + cor_cf*u_perp + sgrad_cf
rhs_con_vector = pinvJ_cf ⋅ rhs_vector # exact contravariant component

# check geostropohic balance
geo_balance = cor_cf*u_perp + sgrad_cf
geo_balance_con = pinvJ_cf⋅ geo_balance
e_geo_balance = sum(∫( geo_balance_con )dΩ)


  # weak forms
  detg_cf = panelwise_cellfield(detg,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)

  u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_contra_h = interpolate(u_contra_cf,U)
  h_h = interpolate(h_cf,P)

  Aperp = [0 0 0
          0 0 -1
          0 1 0]
  Rperp = TensorValue(Aperp)
  Rperp_cf = CellField(Rperp,Ω_panel)

  #### Velocity
  biformU(u,v) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( ( cor_cf*( (Rperp_cf⋅ u)⋅v))*detg_cf )dΩ
  liformU(v) = ∫( rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ + ∫( h_h*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
  op = AffineFEOperator(biformU,liformU,U,V)
  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  uh = FEFunction(U,x)

  uh_proj = covarient_basis_cf ⋅ uh
  e_u = l2( (u_proj_cf - uh_proj),meas_cf,dΩ) # error in physical velocity u





  biformP(p,q) = ∫( (p*q)*meas_cf )dΩ
  liformP(q) =  ∫( (rhs_scalar*q)*meas_cf )dΩ - ∫( q*(u_contra_h⋅grad_meas_cf + meas_cf*(∇⋅u_contra_h) )  )dΩ
  op = AffineFEOperator(biformP,liformP,P,Q)

  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  ph = FEFunction(P,x)

  e_p = l2((h_cf - ph),meas_cf,dΩ) # error in depth


cell_geo_map = geo_map_func(Ω_panel)
panel_cfs = [h_cf, u_proj_cf, ph, uh_proj, h_cf-ph, u_proj_cf-uh_proj ]
labels = ["p","u_proj", "ph", "uh_proj", "ep","eu"]


cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/linear_wave_sol",cellfields=cellfields,append=false,geo_map=cell_geo_map)




return e_u,e_p,false
end
