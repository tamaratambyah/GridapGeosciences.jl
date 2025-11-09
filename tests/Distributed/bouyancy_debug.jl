using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using Gridap.Algebra

using DrWatson
dir = datadir("DistributedLinearisedBoussinesq")
!isdir(dir) && mkdir(dir)

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

include("../convergence_tools.jl")
include("missing_overloads.jl")

function b0(xyz)
  x,y,z = xyz
  k = sqrt(x^2 + y^2 + z^2) - _R

  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr

  θc = 2*π/3
  ϕc = 0.0

  r = _R*acos( sin(ϕc)*sin(ϕ) + cos(ϕc)*cos(ϕ)*cos(θ-θc)    )
  s = _d^2/(_d^2 + r^2)
  b = s*sin(2*π*k/_Lz)
  b

end

function u0(xyz)
  x,y,z = xyz
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr

  # u = 0.04*cos(ϕ) #
  # v = 0.0#
  u = -1.0*_u_0*y
  v = _u_0*x

  VectorValue(u,v,0.0)
end

function un(xyz)
  u = u0(xyz)
  n = normal_vec(xyz)
  u⋅n
end

models = get_3D_octree_refined_models(ranks,4)

# models = get_3D_octree_vertical_refined_models(ranks,5,3)
panel_model = models[2]
p_fe = 1
ls = LUSolver()
return_vtk = true


b = panel_to_cartesian(b0)
vX = panel_to_cartesian(tangent_vec(u0))
_un = panel_to_cartesian(un)
function bouyancy_debug(panel_model,p_fe::Int,dir::String,
    b::Function,vX::Function,_un::Function,
    ls=LUSolver(),return_vtk=false)

    das = FullyAssembledRows()
    # das = SubAssembledRows()

    Dc = num_cell_dims(panel_model)
    println("Dc = $Dc")
    println(num_cells(panel_model))

    lvl_h = nref(nc_horizontal(panel_model))
    lvl_v = nref(nc_vertical(panel_model))
    i_am_main(ranks) && println("nref_h = $lvl_h; nref_v = $lvl_v; p_fe = $p_fe")

    panel_ids = get_panel_ids(panel_model)
    Ω_panel = Triangulation(das,panel_model)
    dΩ = Measure(Ω_panel,4*(p_fe+1))

    _Ω_panel = Triangulation(das,panel_model)
    _dΩ = Measure(_Ω_panel,8*(p_fe+1))

    b_cf = panelwise_cellfield(b,Ω_panel,panel_ids)

    tags = ["bottom_boundary",  "top_boundary"]

    R = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe);
      conformity=:L2,dirichlet_tags=tags)
    B = TrialFESpace(R,b_cf)

    # weak forms
    b_h = interpolate(b_cf,B)

    meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
    metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
    u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
    un_cf = panelwise_cellfield(_un,Ω_panel,panel_ids)

    # rhs_bouyancy = b_cf
    rhs_bouyancy = b_cf + _N^2*un_cf

    n_cf = CellField(VectorValue(1,0,0),Ω_panel)

    assem = SparseMatrixAssembler(B,R,das)

    biformB(b,r) = ∫( (b*r)*meas_cf )dΩ
    liformB(r) =  ∫( (rhs_bouyancy*r)*meas_cf )dΩ - ∫( _N^2*( r*( n_cf⋅(metric_cf⋅u_contra_cf))*meas_cf)   )dΩ

    op = AffineFEOperator(biformB,liformB,B,R,assem)


    A = get_matrix(op)
    b_vec = get_vector(op)
    ns = numerical_setup(symbolic_setup(ls,A),A)
    x = allocate_in_domain(A); fill!(x,0.0)
    solve!(x,ns,b_vec)
    bh = FEFunction(B,x)

    l2(e,meas_cf,dΩ) = sum(∫( (e⋅e)*meas_cf )dΩ)
    e_b = l2((b_h - bh),meas_cf,_dΩ) # error in bouyancy


  if return_vtk
    cell_geo_map = geo_map_func(get_panel_ids(Ω_panel))
    panel_cfs = [b_cf, b_h, bh,  b_cf-bh]
    labels = [ "b_cf", "b_h_int", "bh",  "eb"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nrefh$(lvl_h)_nrefv$(lvl_v)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

    return e_b,false,false
end



# errors,ns,dxs,slopes = h_convergence_test_vertical(models,bouyancy_debug,p_fe,dir,
#     b,ls,true)
errors,ns,dxs,slopes = h_convergence_test(models,bouyancy_debug,p_fe,dir,
    b,vX,_un,ls,true)
plot()
plot_convergence(errors,ns,dxs,slopes;leginf=["b"],colors=[palette(:tab10)[p_fe] ] )
plot!(show=true)
savefig(dir*"/convergence_linear_boussineseq_3D")
