using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using Gridap.ReferenceFEs, Gridap.Polynomials, Gridap.CellData

using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Test
using GridapPETSc
using LinearAlgebra

include("../convergence_tools.jl")

inv_jacobian(p) = x -> inv(forward_jacobian_3D(p)(x))
contra_v_3D(vecX::Function,p::Int) = x -> inv_jacobian(p)(x) ⋅ vecX(p)(x)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)

transpose_jacobian(p) = x -> transpose(forward_jacobian_3D(p)(x))
inv_tranpose_jacobian(p) = x -> inv(transpose_jacobian(p)(x))

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

dir = datadir("HodgeLaplacianConvergence")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

n_ref_lvls = 3
ps = [0,1]
radius = 1.0
thickness = 0.19
models  = get_3D_octree_refined_models(ranks,n_ref_lvls,radius,thickness)
ls = LUSolver()

panel_model = models[2]
p_fe = ps[2]

degree = 4*(p_fe + 1)
if p_fe == 0
  degree = 8
end
@check degree > 0 "Zero quad!!"

## finite element solver
panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,2*degree)
Ω_error = Triangulation(panel_model)
dΩ_error = Measure(Ω_error,4*degree)

tags = ["top_boundary", "bottom_boundary"]
Γ = BoundaryTriangulation(panel_model,tags=tags)
dΓ = Measure(Γ,2*degree)
nΓ = get_normal_vector(Γ)

## metric information
inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω_panel,panel_ids)



function uX(p)
  function _u(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    VectorValue(-xyz[2],xyz[1],0.0)

    # r = sqrt(xyz[1]^2 + xyz[2]^2 + xyz[3]^2)
    # f = 2.0*xyz[3]/r
    # n = normal_vec(xyz)
    # f*n
  end
end


function ucov(p)
  function _u(γαβ)
    u = uX(p)(γαβ)
    J = transpose_jacobian(p)(γαβ)
    J⋅u
  end
end

function unX(p)
  function _u(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    u = uX(p)(γαβ)
    n = normal_vec(xyz)
    u⋅n
  end
end

curlu(p,x) = curl(ucov(p))(x)
curlu(p) = x -> curlu(p,x)
_curlu(p) = curl(ucov(p))


wcov(p,x) = 1/sqrtg(p,x)*metric(p,x)⋅curlu(p,x)
wcov(p) = x -> wcov(p,x)
curlw(p,x) = curl(wcov(p))(x)
curlw(p) = x -> curlw(p,x)
curlw(1)(Point(1,1,1))

curlw_cov(p) = x -> 1/sqrtg(p,x)*metric(p,x)⋅curlw(p,x)
curlw_cov(1)(Point(1,1,1))

_area_meas(p) = x->  forward_jacobian_3D(p,x) ⋅ (inv_metric(p,x) ⋅ VectorValue(1,0,0))
area_meas(p) = x-> norm(_area_meas(p)(x))
area_meas(1)(Point(1,1,1))

wcrossk_cov(p) = x -> 1/sqrtg(p,x) * metric(p,x)⋅(wcov(p,x) × (VectorValue(1,0,0)/area_meas(p)(x)) )
wcrossk_cov(1)(Point(1,1,1))





_sdiv_u(p) = x -> sqrtg(p)(x)* contra_v_3D(uX,p)(x)
sdiv_u(p) = x -> 1/sqrtg(p)(x) * ( divergence(_sdiv_u(p) )(x) )
graddiv_cov(p) = gradient(sdiv_u(p))
rhs(p) = x-> curlw_cov(p)(x) - graddiv_cov(p)(x)

rhs_cov_cf = panelwise_cellfield(rhs,Ω_panel,panel_ids)

u_cov_cf = panelwise_cellfield(ucov,Ω_panel,panel_ids)
ccurlu_cov_cf = panelwise_cellfield(curlw_cov,Ω_panel,panel_ids)
un_cf = panelwise_cellfield(unX,Ω_panel,panel_ids)
curlu_cross = panelwise_cellfield(wcrossk_cov,Ω_panel,panel_ids)
curlu_cf = panelwise_cellfield(curlu,Ω_panel,panel_ids)

sdiv_cf =  panelwise_cellfield(surfdiv(contra_v_3D(uX)),Ω_panel,panel_ids)
sigma_cf = -sdiv_cf


cellfields = ["curlu"=>ccurlu_cov_cf,
              "u"=>covariant_basis_cf ⋅ (inv_metric_cf⋅u_cov_cf),
              "un"=>un_cf,
              "curlu_cross"=>covariant_basis_cf ⋅ (inv_metric_cf⋅curlu_cross),
              "sigma"=>-sdiv_cf,
              "rhs"=>covariant_basis_cf ⋅ (inv_metric_cf⋅ rhs_cov_cf)
               ]
writevtk(Ω_panel,dir*"/sol",
        cellfields=cellfields,
        append=false,geo_map= geo_map_func(Ω_panel))





## FE spaces
T = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
S = TrialFESpace(T)

## Sigma
biform_s(s,t) =  ∫( (s*t)*meas_cf  )dΩ
boundary(t) = (
              ∫( ( t*(u_cov_cf⋅(inv_metric_cf⋅nΓ)) )*(meas_cf)  )dΓ
                    )
liform_s(t) = ( ∫( ∇(t)⋅(inv_metric_cf⋅u_cov_cf)*meas_cf  )dΩ
                - boundary(t)
                )

# boundary_vec = assemble_vector(boundary,T)
# partition(boundary_vec).item
# norm(boundary_vec)


op = AffineFEOperator(biform_s,liform_s,S,T)
sh = solve(ls,op)

# A = get_matrix(op)
# eigvals(Array(partition(A).item))


_e = sigma_cf - sh
el2_s = sqrt(sum(∫( (_e*_e)*meas_cf  )dΩ_error))

cellfields = ["sh"=>sh, "s"=>sigma_cf, "e"=>sh-sigma_cf
               ]
writevtk(Ω_panel,dir*"/sol",
        cellfields=cellfields,
        append=false,geo_map= geo_map_func(Ω_panel))





## u
R = TestFESpace(Ω_panel, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl)
H = TrialFESpace(R)

biform_u(u,v) = ( ∫( u⋅(inv_metric_cf⋅v)*meas_cf  )dΩ
                +   ∫( curl(u)⋅(metric_cf⋅curl(v))*(1/meas_cf) )dΩ
                )
liform_u(v) = ( ∫( u_cov_cf⋅(inv_metric_cf⋅v)*meas_cf )dΩ
              + ∫( rhs_cov_cf⋅(inv_metric_cf⋅v)*meas_cf  )dΩ
              + ∫( -1.0*gradient(sigma_cf)⋅(inv_metric_cf⋅v)*meas_cf )dΩ
              + ∫( v⋅( ( metric_cf⋅curlu_cf )×nΓ    )*(1/meas_cf)     )dΓ
                )
op = AffineFEOperator(biform_u,liform_u,H,R)
uh = solve(ls,op)

# A = get_matrix(op)
# evl = real.(eigvals(Array(partition(A).item)))

_e = (inv_metric_cf⋅uh) - (inv_metric_cf⋅u_cov_cf)
el2_u = sqrt(sum(∫( (_e⋅(metric_cf ⋅_e))*meas_cf  )dΩ_error))



cellfields = ["u"=>covariant_basis_cf ⋅ (inv_metric_cf⋅u_cov_cf),
"uh"=>covariant_basis_cf ⋅ (inv_metric_cf⋅uh),
"eu"=>covariant_basis_cf ⋅ (inv_metric_cf⋅uh)-covariant_basis_cf ⋅ (inv_metric_cf⋅u_cov_cf),
               ]
writevtk(Ω_panel,dir*"/sol",
        cellfields=cellfields,
        append=false,geo_map= geo_map_func(Ω_panel))


### Multifield
X = MultiFieldFESpace([S,H])
Y = MultiFieldFESpace([T,R])

biform_x((s,u),(t,v)) = (
                ∫( (s*t)*meas_cf  )dΩ
              - ∫( ∇(t)⋅(inv_metric_cf⋅u)*meas_cf  )dΩ
              + ∫( curl(u)⋅(metric_cf⋅curl(v))*(1/meas_cf) )dΩ
              + ∫( gradient(s)⋅(inv_metric_cf⋅v)*meas_cf )dΩ
                )
liform_x((t,v)) = (
               ∫( rhs_cov_cf⋅(inv_metric_cf⋅v)*meas_cf  )dΩ
              + ∫( v⋅( ( metric_cf⋅curlu_cf )×nΓ    )*(1/meas_cf)     )dΓ
              - ∫(( t*(u_cov_cf⋅(inv_metric_cf⋅nΓ)) )*(meas_cf)  )dΓ
                )


op = AffineFEOperator(biform_x,liform_x,X,Y)
sh,uh = solve(ls,op)


_e = sigma_cf - sh
el2_s = sqrt(sum(∫( (_e*_e)*meas_cf  )dΩ_error))

_e = (inv_metric_cf⋅uh) - (inv_metric_cf⋅u_cov_cf)
el2_u = sqrt(sum(∫( (_e⋅(metric_cf ⋅_e))*meas_cf  )dΩ_error))


cellfields =  ["u"=>covariant_basis_cf ⋅ (inv_metric_cf⋅u_cov_cf),
"uh"=>covariant_basis_cf ⋅ (inv_metric_cf⋅uh),
"eu"=>covariant_basis_cf ⋅ (inv_metric_cf⋅uh)-covariant_basis_cf ⋅ (inv_metric_cf⋅u_cov_cf),
"sh"=>sh, "s"=>sigma_cf, "e"=>sh-sigma_cf
               ]
writevtk(Ω_panel,dir*"/sol",
        cellfields=cellfields,
        append=false,geo_map= geo_map_func(Ω_panel))
