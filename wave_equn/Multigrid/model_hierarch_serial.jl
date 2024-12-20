using GridapGeosciences
using Gridap
using DrWatson
using FillArrays
using GridapSolvers
using Gridap.Helpers: @check
import GridapSolvers.MultilevelTools: ModelHierarchyLevel, HierarchicalArray
import Gridap.Adaptivity: EdgeBasedRefinement, RedGreenRefinement
import Gridap.Geometry: UnstructuredDiscreteModel




function Gridap.Adaptivity.refine(method::EdgeBasedRefinement,model::UnstructuredDiscreteModel{Dc,Dp};cells_to_refine=nothing) where {Dc,Dp}
  println("refined method")
  function _cell_vector_to_dof_vector!(dof_vector,cell_node_ids, cell_vector)
    cache_cell_node_ids = array_cache(cell_node_ids)
    cache_cell_vector   = array_cache(cell_vector)
    for k=1:length(cell_node_ids)
       current_node_ids = getindex!(cache_cell_node_ids,cell_node_ids,k)
       current_values   = getindex!(cache_cell_vector,cell_vector,k)
       for (i,id) in enumerate(current_node_ids)
        dof_vector[current_node_ids[i]]=current_values[i]
       end
    end
  end


  ctopo = Gridap.Geometry.get_grid_topology(model)
  coarse_labels = get_face_labeling(model)
  # Create new model
  rrules, faces_list = Gridap.Adaptivity.setup_edge_based_rrules(method, model.grid_topology,cells_to_refine)
  topo   = Gridap.Adaptivity.refine_edge_based_topology(ctopo,rrules,faces_list)

  reffes = map(p->Gridap.FESpaces.LagrangianRefFE(Float64,p,1),Gridap.Geometry.get_polytopes(topo))
  grid   = Gridap.Geometry.UnstructuredGrid(
    Gridap.Geometry.get_vertex_coordinates(topo),
    Gridap.Geometry.get_faces(topo,Dc,0),reffes,
    Gridap.Geometry.get_cell_type(topo),
    Gridap.Geometry.OrientationStyle(topo)
  )
  glue = Gridap.Adaptivity.blocked_refinement_glue(rrules)
  labels = Gridap.Adaptivity.refine_face_labeling(coarse_labels,glue,model.grid_topology,topo)
  ref_model = Gridap.Geometry.UnstructuredDiscreteModel(grid,topo,labels)

  # adapted_ref_model = Gridap.Adaptivity.AdaptedDiscreteModel(ref_model,model,glue)
  # writevtk(adapted_ref_model,joinpath(datadir("models"),"ref_model"), append=false)


  cube_surface_trian = Triangulation(ref_model)
  order = 1
  radius = 1

  vector_reffe=ReferenceFE(lagrangian,VectorValue{3,Float64},order)
  V = FESpace(cube_surface_trian,vector_reffe; conformity=:H1)
  vh = interpolate(GridapGeosciences.MapCubeToSphere(radius),V)

  scalar_reffe=ReferenceFE(QUAD,lagrangian,Float64,order)
  xref=Gridap.ReferenceFEs.get_node_coordinates(scalar_reffe)

  xrefₖ=Fill(xref,num_cells(ref_model))

  vhx=lazy_map(evaluate,Gridap.CellData.get_data(vh),xrefₖ)

  V = FESpace(cube_surface_trian,scalar_reffe; conformity=:H1)
  node_coordinates = Vector{Point{3,Float64}}(undef,num_free_dofs(V))
  cell_node_ids    = get_cell_dof_ids(V)
  _cell_vector_to_dof_vector!(node_coordinates,cell_node_ids,vhx)
  cell_types  = collect(Fill(1,num_cells(cube_surface_trian)))
  cell_reffes = [scalar_reffe]

  cube_surface_grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                                Gridap.Arrays.Table(cell_node_ids),
                                                  cell_reffes,
                                                  cell_types) ### NOT ORIENTED!
  cube_model = Gridap.Geometry.UnstructuredDiscreteModel(cube_surface_grid,topo,labels)

  # adapted_model = Gridap.Adaptivity.AdaptedDiscreteModel(cube_model,model,glue)
  # writevtk(adapted_model,joinpath(datadir("models"),"cube_model"), append=false)

  return Gridap.Adaptivity.AdaptedDiscreteModel(cube_model,model,glue)

end

# coarsest model
_model = CubedSphereDiscreteModel(2)
model = _model.cubed_sphere_linear_model
model2 = Gridap.Adaptivity.refine(RedGreenRefinement(),model;cells_to_refine=nothing)#Gridap.Adaptivity.refine(model, refinement_method = "red_green")
model3 = Gridap.Adaptivity.refine(model2, refinement_method = "red_green")

models = [model3,model2]

# writevtk(model,joinpath(datadir("models"),"cubed_sphere_serial_linear"), append=false)
# writevtk(model2,joinpath(datadir("models"),"cubed_sphere_serial_refine"), append=false)

mh = ModelHierarchy(models)

#### test grid transfer
qorder = 6
order = 0
sol(x) = sum(x)
# Triangulations
trian = Triangulation(model2)
ctrian = Triangulation(model)

glue = Gridap.Adaptivity.get_adaptivity_glue(model2)
rrules = Gridap.Adaptivity.get_old_cell_refinement_rules(glue)

# Measures
dΩ_f  = Measure(trian,qorder)
dΩ_c  = Measure(ctrian,qorder)
dΩ_cf = Measure(ctrian,trian,qorder)

cell_quad = Gridap.CellData.get_cell_quadrature(dΩ_cf)
dΩ_cf_bis = Measure(ctrian,trian,cell_quad)

# FESpaces
reffe = ReferenceFE(lagrangian,Float64,order)
V_f = TestFESpace(model2,reffe;dirichlet_tags="boundary")
U_f = TrialFESpace(V_f,sol)
V_c = TestFESpace(model,reffe;dirichlet_tags="boundary")
U_c = TrialFESpace(V_c,sol)
# V_c_fast = TestFESpace(model,rrules,reffe;dirichlet_tags="boundary")
# U_c_fast = TrialFESpace(V_c_fast,sol)

# CellField: Coarse -> Fine
cf_c_phy = CellField(sol,ctrian)
cf_c_ref = Gridap.CellData.change_domain(cf_c_phy,PhysicalDomain(),ReferenceDomain())
cf_f_ref_ref = Gridap.CellData.change_domain(cf_c_ref, trian, ReferenceDomain())
cf_f_ref_phy = Gridap.CellData.change_domain(cf_c_ref, trian, PhysicalDomain())
cf_f_phy_ref = Gridap.CellData.change_domain(cf_c_phy, trian, ReferenceDomain())
cf_f_phy_phy = Gridap.CellData.change_domain(cf_c_phy, trian, PhysicalDomain())




function get_jacobi_smoothers(mh)
  nlevs = num_levels(mh)
  smoothers = Fill(RichardsonSmoother(JacobiLinearSolver(),10,2.0/3.0),nlevs-1)
  level_parts = view(get_level_parts(mh),1:nlevs-1)
  return HierarchicalArray(smoothers,level_parts)
end

function get_patch_smoothers(mh,tests,biform,qdegree,ω=0.2,niter=10,local_solver=LUSolver(),is_nonlinear=false,weighted=false)
  patch_decompositions = PatchDecomposition(mh)
  patch_spaces = PatchFESpace(tests,patch_decompositions)
  nlevs = num_levels(mh)
  smoothers = map(view(tests,1:nlevs-1),patch_decompositions,patch_spaces) do tests, PD, Ph
    Vh = get_fe_space(tests)
    Ω  = Triangulation(PD)
    dΩ = Measure(Ω,qdegree)
    ap = (u,v) -> biform(u,v,dΩ)
    patch_smoother = GridapSolvers.PatchBasedSmoothers.PatchBasedLinearSolver(ap,Ph,Vh;local_solver,is_nonlinear,weighted)
    return RichardsonSmoother(patch_smoother,niter,ω) # 10 iterations, ω = 0.2
  end
  return smoothers
end

function get_bilinear_form(mh_lev,biform,qdegree)
  model = get_model(mh_lev)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,qdegree)
  return (u,v) -> biform(u,v,dΩ)
end

u_exact(x) = VectorValue(x[1]*(1-x[1]), x[2]*(1-x[2]),0.0)
p_exact(x) = 1.0 + 0.001*x[1]*(1-x[1]) + 0.001*x[2]*(1-x[2])
f_u(x) =  u_exact(x) + ∇(p_exact)(x)
f_p(x) = p_exact(x) + (∇⋅u_exact)(x)


p = 1
degree = 6



model = get_model(mh,1)
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
l2(e,dΩ) = sum(∫(e⊙e)dΩ)

reffe_u = ReferenceFE(raviart_thomas,Float64,p)
Vs = TestFESpace(mh,reffe_u,conformity=:Hdiv)
Us = TrialFESpace(Vs)

reffe_p = ReferenceFE(lagrangian,Float64,p)
Qs = TestFESpace(mh,reffe_p;conformity=:L2)
Ps = TrialFESpace(Qs)

tests = MultiFieldFESpace([Us,Ps])
trials = MultiFieldFESpace([Vs,Qs])


biform((u,p),(v,q),dΩ) = ( ∫( u⋅v - (∇⋅v)*p )dΩ
                         + ∫( p*q + (∇⋅u)*q )dΩ )
liform((v,q),dΩ) = ∫( f_u⋅v  + f_p⋅q )dΩ

a(u,v) = biform(u,v,dΩ)
l(v) = liform(v,dΩ)
op = AffineFEOperator(a,l,get_fe_space(trials,1),get_fe_space(tests,1))
A = get_matrix(op)
b = get_vector(op)

biforms = map(mhl -> get_bilinear_form(mhl,biform,degree),mh)

smoothers = get_patch_smoothers(mh,tests,biform,degree)
# smoothers = get_jacobi_smoothers(mh)

restrictions, prolongations = setup_transfer_operators(
  trials, degree; mode=:residual, solver=GridapSolvers.LinearSolvers.IS_ConjugateGradientSolver(;reltol=1.e-6)
)


gmg = GMGLinearSolver(
  mh,trials,tests,biforms,
  prolongations,restrictions,
  pre_smoothers=smoothers,
  post_smoothers=smoothers,
  coarsest_solver=LUSolver(),
  maxiter=3,mode=:preconditioner,
  verbose=true
)
gmg.log.depth = 2


# smatrices, A, b = compute_hierarchy_matrices(trials,tests,biform,liform,degree)

solver = GMRESSolver(20;Pl=gmg,maxiter=5,atol=1e-14,rtol=1.e-14,verbose=true)
ns = numerical_setup(symbolic_setup(solver,A),A)

# Solve
x = Gridap.Algebra.allocate_in_domain(A); fill!(x,0.0)
solve!(x,ns,b)


# Error
Uh = get_fe_space(trials,1)
uh, ph = FEFunction(Uh,x)

l2(u_exact-uh,dΩ)
l2(p_exact-ph,dΩ)

writevtk(Ω,datadir("wave"),cellfields = ["uh"=>uh, "ph"=>ph,
                "u"=>CellField(u_exact,Ω), "p"=>CellField(p_exact,Ω)],append=false)
