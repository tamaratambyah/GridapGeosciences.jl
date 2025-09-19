using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

include("AdvectionSUPG.jl")
include("AdvectionDGUpwinding.jl")

#### Replicate test in Section 5.4 of Rognes2013 paper
vecX(XYZ) = VectorValue(-XYZ[2],XYZ[1],0.0)
u0(XYZ) = exp(-(XYZ[2]^2 + XYZ[3]^2)  )
u0vecX(XYZ) = u0(XYZ)*vecX(XYZ)

vX = panel_to_cartesian(tangent_vec(vecX))
u = panel_to_cartesian(u0)
uvX = panel_to_cartesian(u0vecX)

n_ref_lvls = 4
CFL = 0.1


## SUPG
advection_supg_convergence_test(n_ref_lvls,u,vX,CFL,true)
transient_advection_supg_convergence_test(n_ref_lvls,u,vX,CFL,false)

## DG
advection_dg_convergence_test(n_ref_lvls,u,vX,uvX,true)
transient_advection_dg_convergence_test(n_ref_lvls,u,vX,CFL,false)
