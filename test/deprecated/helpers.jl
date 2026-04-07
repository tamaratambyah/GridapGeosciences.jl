import GridapGeosciences.Helpers: RADIUS
import GridapGeosciences.Helpers: THICKNESS

l2(e,dΩ) = sum(∫( e⋅e )dΩ)
l2(e,meas,dΩ) = sum(∫( (e⋅e)*meas )dΩ)
mass_conservation(p,meas,dΩ) = sum(∫( (p)*meas )dΩ)
# dx(nc) = sqrt( 4*π*RADIUS^2 / (6*sqrt(nc)^2) )

nc(model::ParametricOctreeDistributedDiscreteModel) = nc(model.parametric_dmodel)
dx(model::ParametricOctreeDistributedDiscreteModel) = dx(model.parametric_dmodel)
