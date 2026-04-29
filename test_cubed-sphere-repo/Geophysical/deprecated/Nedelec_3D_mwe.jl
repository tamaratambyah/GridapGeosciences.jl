using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using MPI
using PartitionedArrays
using MPIPreferences
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using GridapGeosciences
using GridapPETSc


# include("CurlConformingFESpacesFixes.jl")


# Initial fluid depth
function h_3D(p)
  function _h(γαβ)
    ζ = 0.0
    xyz = forward_map_3D(p)(γαβ)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    h  = -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ)
    _H_0 - (_ω*_u0 + 0.5*_u0*_u0)*h*h/_g
  end
end

function f_3D(p)
  function _f(γαβ)
    ζ = 0.0
    xyz = forward_map_3D(p)(γαβ)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    2.0*_ω*( -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ) )
  end
end

function f_vec_3D(p)
  function _f(γαβ)
    f = f_3D(p)(γαβ)

    xyz = forward_map_3D(p)(γαβ)
    n = normal_vec(xyz)
    f*n
  end
end



function η_3D(p)
  function _η(γαβ)
    ζ = 0.0
    xyz = forward_map_3D(p)(γαβ)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    η = -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ)
    ( 4*π/_T + 2*_ω  )* η
  end
end


function η_vec_3D(p)
  function _η(γαβ)
    η = η_3D(p)(γαβ)

    xyz = forward_map_3D(p)(γαβ)
    n = normal_vec(xyz)
    η*n
  end
end


function u_vec_3D(p)
  function _u(γαβ)
    ζ = 0.0
    xyz = forward_map_3D(p)(γαβ)
    # θϕr   = xyz2θϕr(xyz)
    # θ,ϕ,r = θϕr
    # u     = _u0*(cos(ϕ)*cos(ζ) + cos(θ)*sin(ϕ)*sin(ζ))
    # v     = - _u0*sin(θ)*sin(ζ)
    # _spherical_to_cartesian_matrix(θϕr)⋅VectorValue(u,v,0)

    #### from Rognes2013 paper
    r = sqrt(xyz[1]^2 + xyz[2]^2 + xyz[3]^2)
    u = -_u0*xyz[2]/r
    v = _u0*xyz[1]/r
    w = 0.0
    VectorValue(u,v,w)
  end
end

THICKNESS = 3e-4

#### Parameters
a_e = 6.37e6 # m
g = 9.8 # m/2
ω = 7.29e-5 #s^-1
T = 12*24*3600 #s
H_0 = 5960.0 #2.94e4/g #m
u_0 = 20 #2*π*a_e/T #m/s
b_0 = 2000.0 #m
TF = 20*(3600*24)

## thickness computation
b_0 = 2000 #m
n_ref = 4
n = 6*2^n_ref*2^n_ref
dxx = sqrt(4*π*a_e^2/n)
dz = b_0

s = dxx/dz
_dx = sqrt(4*π/n)
t = _dx/s


L = a_e
_τ = 1/ω
LV = b_0/THICKNESS

_a = a_e/L
_g = g*_τ^2/L
_ω = ω*_τ
_H_0 = H_0/L
_T = T/_τ
_u0 = u_0/L*_τ
_b0 = b_0/LV
_tF = TF/_τ


inv_jacobian(p) = x -> inv(forward_jacobian(p)(x))
contra_v_3D(vecX::Function,p::Int) = x -> inv_jacobian(p)(x) ⋅ vecX(p)(x)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)

transpose_jacobian(p) = x -> transpose(forward_jacobian(p)(x))
inv_tranpose_jacobian(p) = x -> inv(transpose_jacobian(p)(x))
contravariant_basis_3D(p) = x -> inv_tranpose_jacobian(p)(x)

covar_v_3D(vecX::Function,p::Int) = x -> transpose_jacobian(p)(x) ⋅ vecX(p)(x)
covar_v_3D(vecX::Function) = p -> covar_v_3D(vecX,p)


MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
ranks = distribute_with_mpi(LinearIndices((nprocs,)))

ls = LUSolver()
p_fe = 1
n_ref_lvls = 3

o3model = GridapGeosciences.Distributed.CubedSphere3DParametricOctreeDistributedDiscreteModel(ranks;
    num_horizontal_uniform_refinements=n_ref_lvls,
    num_vertical_uniform_refinements=0)
panel_model = o3model.parametric_dmodel



## finite element solver
panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,4*(p_fe+1))

tags = ["top_boundary", "bottom_boundary"]
Γ = BoundaryTriangulation(panel_model,tags=tags)
dΓ = Measure(Γ,4*(p_fe+1))
nΓ = get_normal_vector(Γ)

R = TestFESpace(Ω_panel, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
H = TrialFESpace(R,VectorValue(0.0,0.0,0.0))

Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)

V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

X_prog = MultiFieldFESpace([U,P]) # u, p
Y_prog = MultiFieldFESpace([V,Q]) # u, p

## initial conditions
u_contra_cf = ParametricCellField(contra_v_3D(u_vec_3D),Ω_panel,panel_ids)
u_contra_h = interpolate(u_contra_cf,U)

h_cf = ParametricCellField(h_3D,Ω_panel,panel_ids)
h_h = interpolate(h_cf,P)

xh0 = interpolate_everywhere([u_contra_h,h_h],X_prog)

inv_metric_cf = ParametricCellField(inv_metric,Ω_panel,panel_ids)
metric_cf = ParametricCellField(metric,Ω_panel,panel_ids)
meas_cf = ParametricCellField(sqrtg,Ω_panel,panel_ids)
covariant_basis_cf = ParametricCellField(covariant_basis,Ω_panel,panel_ids)
jac_cf = ParametricCellField(forward_jacobian,Ω_panel,panel_ids)
area_meas_cf = Operation(norm)(jac_cf⋅(inv_metric_cf ⋅nΓ) )

gravity = _g
f_cov_cf = ParametricCellField(covar_v_3D(f_vec_3D),Ω_panel,panel_ids)


uh,ph = xh0
t0 = 0.0

### single field
biformq(q,w) = ∫( ph*(q⋅(inv_metric_cf⋅w))*meas_cf )dΩ
liformq(w) = (
              ∫( uh⋅( metric_cf⋅ curl(w) )  )dΩ
            - ∫( (( w × (metric_cf⋅ uh) )⋅nΓ)*area_meas_cf   )dΓ
            + ∫( (f_cov_cf⋅(inv_metric_cf ⋅ w))*meas_cf )dΩ
              )
op = AffineFEOperator(biformq,liformq,H,R)
qh = solve(ls,op)


vort = qh*ph
vortf = vort - f_cov_cf

vort_ambient = covariant_basis_cf ⋅ (inv_metric_cf ⋅ vort )
vortf_ambient = covariant_basis_cf ⋅ (inv_metric_cf ⋅ vortf )
q_ambient = covariant_basis_cf ⋅ (inv_metric_cf ⋅ qh )

cellfields = ["q"=>q_ambient, "vort"=>vort_ambient, "vortf"=>vortf_ambient ]
dir = datadir("Nedelec_3D")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)
writevtk(Ω_panel,dir*"/IC_nrproc$(nprocs)",cellfields=cellfields,append=false,geo_map= latlon_geo_map_func(Ω_panel))
