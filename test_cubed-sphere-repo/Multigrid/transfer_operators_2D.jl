using DrWatson
using PartitionedArrays
using Gridap, Gridap.Algebra
using Gridap.Helpers, Gridap.Adaptivity
using Test

using GridapSolvers
using GridapSolvers.MultilevelTools
using GridapSolvers.LinearSolvers
using GridapGeosciences

using MPI, PartitionedArrays, GridapP4est, GridapDistributed
using Plots

include("../Laplace/analytic_funcs.jl")
include("../convergence_tools.jl")

# models = get_refined_models(3,true)
# mh = ModelHierarchy(models)

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

n_ref_lvls = 3
radius = 1
model0 = ParametricOctreeDistributedDiscreteModel(ranks, radius; num_initial_uniform_refinements=0)
mh = ModelHierarchy(model0,n_ref_lvls)

f = panel_to_cartesian(fθϕ)

# FE Spaces
p_fe  = 1
reffe  = ReferenceFE(lagrangian,Float64,p_fe)
tests  = TestFESpace(mh,reffe;conformity=:H1)
trials = TrialFESpace(tests)
solver = LUSolver()

nlevs = num_levels(mh)

degress = [i*(2*p_fe+1) for i in 1:4]

errors_c2f_projection = [Vector{Float64}(undef,nlevs-1) for i in 1:length(degress)]
errors_f2c_projection = [Vector{Float64}(undef,nlevs-1) for i in 1:length(degress)]

errors_c2f_interpolation = [Vector{Float64}(undef,nlevs-1) for i in 1:length(degress)]
errors_f2c_interpolation = [Vector{Float64}(undef,nlevs-1) for i in 1:length(degress)]

errors_c2f2c_interpolation = [Vector{Float64}(undef,nlevs-1) for i in 1:length(degress)]
errors_c2f2c_projection = [Vector{Float64}(undef,nlevs-1) for i in 1:length(degress)]


for (i,degree) in enumerate(degress)

  for lev in 1:nlevs-1
    ### Fine mesh
    model_h = get_model_before_redist(mh,lev)
    panel_idsh = get_panel_ids(model_h)
    Vh  = get_fe_space_before_redist(tests,lev)
    Uh  = get_fe_space_before_redist(trials,lev)
    Ωh  = get_triangulation(model_h)
    dΩh = Measure(Ωh,degree)
    dΩh_error = Measure(Ωh,6*degree)
    sol_cf = panelwise_cellfield(f,Ωh,panel_idsh)
    uh  = interpolate(sol_cf,Uh)
    meas_h = panelwise_cellfield(sqrtg,Ωh,panel_idsh)

    ### Coarse mesh
    model_H = get_model(mh,lev+1)
    panel_idsH = get_panel_ids(model_H)
    VH  = get_fe_space(tests,lev+1)
    UH  = get_fe_space(trials,lev+1)
    ΩH  = get_triangulation(model_H)
    dΩH = Measure(ΩH,degree)
    dΩH_error = Measure(ΩH,6*degree)
    sol_cf = panelwise_cellfield(f,ΩH,panel_idsH)
    uH  = interpolate(sol_cf,UH)
    meas_H = panelwise_cellfield(sqrtg,ΩH,panel_idsH)

    dΩhH = Measure(ΩH,Ωh,degree)

    ############################################################################
    ##### Coarse FEFunction -> Fine FE function
    ############################################################################
    # By interpolation
    uH_i = interpolate(uH,Uh)
    _eh = uH_i-uH
    eh  = sum(∫( (_eh*_eh)*meas_h)*dΩh_error)
    # eh  = sum(∫( (_eh*_eh))*dΩh_error)
    errors_c2f_interpolation[i][lev] = sqrt(eh)

    # By projection
    ah(u,v) = ∫( (v*u)*meas_h  )*dΩh
    lh(v)   = ∫( (v*uH)*meas_h )*dΩh
    oph = AffineFEOperator(ah,lh,Uh,Vh)
    Ah  = get_matrix(oph)
    bh  = get_vector(oph)

    xh = allocate_in_domain(Ah); fill!(xh,0.0)
    ns = numerical_setup(symbolic_setup(solver,Ah),Ah)
    solve!(xh,ns,bh)
    uH_projected = FEFunction(Uh,xh)

    _eh = uH-uH_projected
    # eh  = sum(∫( (_eh*_eh))*dΩh_error)
    eh  = sum(∫( (_eh*_eh)*meas_h)*dΩh_error)
    errors_c2f_projection[i][lev] = sqrt(eh)

    ############################################################################
    ##### Fine FEFunction -> coarse FE function
    ############################################################################
    # By interpolation
    uh_i = interpolate(uh,UH)
    _eH = uh_i-uh
    eH  = sum(∫( (_eH*_eH)*meas_H)*dΩH_error)
    # eH  = sum(∫( (_eH*_eH))*dΩH_error)
    errors_f2c_interpolation[i][lev] = sqrt(eH)

    # By projection
    aH(u,v) = ∫( (v*u)*meas_H )*dΩH
    lH(v)   = ∫( (v*uh)*meas_h )*dΩhH
    opH = AffineFEOperator(aH,lH,UH,VH)
    AH  = get_matrix(opH)
    bH  = get_vector(opH)

    xH = allocate_in_domain(AH); fill!(xH,0.0)
    ns = numerical_setup(symbolic_setup(solver,AH),AH)
    solve!(xH,ns,bH)
    uh_projected = FEFunction(UH,xH)

    _eH = uh-uh_projected
    eH  = sum(∫( (_eH*_eH)*meas_H)dΩH_error)
    # eH  = sum(∫( (_eH*_eH))dΩH_error)
    errors_f2c_projection[i][lev] = sqrt(eH)


    ############################################################################
    ##### Coarse -> Fine FEFunction -> Coarse FE function
    ############################################################################
    # By interpolation
    uh_i = interpolate(uH_i,UH)
    _eH = uh_i-uH_i
    eH  = sum(∫( (_eH*_eH)*meas_H)*dΩH_error)
    # eH  = sum(∫( (_eH*_eH))*dΩH_error)
    errors_c2f2c_interpolation[i][lev] = sqrt(eH)

    # By projection
    aH(u,v) = ∫( (v*u)*meas_H )*dΩH
    lH(v)   = ∫( (v*uH_projected)*meas_h )*dΩhH
    opH = AffineFEOperator(aH,lH,UH,VH)
    AH  = get_matrix(opH)
    bH  = get_vector(opH)

    xH = allocate_in_domain(AH); fill!(xH,0.0)
    ns = numerical_setup(symbolic_setup(solver,AH),AH)
    solve!(xH,ns,bH)
    uh_projected = FEFunction(UH,xH)

    _eH = uH_projected-uh_projected
    eH  = sum(∫( (_eH*_eH)*meas_H)dΩH_error)
    # eH  = sum(∫( (_eH*_eH))dΩH_error)
    errors_c2f2c_projection[i][lev] = sqrt(eH)



  end

end


## the errors go through models fine -> coarse, so xaxis is nlevs-1 -> 1
xs = collect(nlevs-1:-1:1)
plot()
_colors = palette(:tab10)
for (i,degree) in enumerate(degress)
  plot!(xs,errors_c2f_projection[i],
      lw = 3,marker=:circle,ms=6,color=_colors[i],ls=:solid,
      label="c2f: quad = $degree")
  plot!(xs,errors_f2c_projection[i],
      lw = 3,marker=:rect,ms=6,color=_colors[i],ls=:dashdot,
      label= "f2c: quad =$degree")
  plot!(xs,errors_c2f2c_projection[i],
      lw = 3,marker=:star5,ms=6,color=_colors[i],ls=:dashdot,
      label= "c2f2c: quad =$degree")
end
plot!(
  yscale=:log10,
      xlabel="level",ylabel="error",
      title="Projection",
      framestyle = :box,
      show=true)
savefig(plotsdir("omodel_projection_transfer_p$p_fe.png"))





plot()
_colors = palette(:tab10)
for (i,degree) in enumerate(degress)
  plot!(xs,errors_c2f_interpolation[i],
      lw = 3,marker=:circle,ms=6,color=_colors[i],ls=:solid,
      label="c2f: quad = $degree")
  plot!(xs,errors_f2c_interpolation[i],
      lw = 3,marker=:rect,ms=6,color=_colors[i],ls=:dashdot,
      label= "f2c: quad =$degree")
  plot!(xs,errors_c2f2c_interpolation[i],
      lw = 3,marker=:star5,ms=6,color=_colors[i],ls=:dashdot,
      label= "c2f2c: quad =$degree")
end
plot!(
  yscale=:log10,
      xlabel="level",ylabel="error",
      title="Interpolation",
      framestyle = :box,
      show=true)
savefig(plotsdir("omodel_interpolation_transfer_p$p_fe.png"))
