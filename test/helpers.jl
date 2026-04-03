import GridapGeosciences.Helpers: RADIUS
import GridapGeosciences.Helpers: THICKNESS

l2(e,dΩ) = sum(∫( e⋅e )dΩ)
l2(e,meas,dΩ) = sum(∫( (e⋅e)*meas )dΩ)
mass_conservation(p,meas,dΩ) = sum(∫( (p)*meas )dΩ)
dx(nc) = sqrt( 4*π*RADIUS^2 / (6*sqrt(nc)^2) )

nc(model::ParametricOctreeDistributedDiscreteModel) = nc(model.parametric_dmodel)
dx(model::ParametricOctreeDistributedDiscreteModel) = dx(model.parametric_dmodel)



function dx(model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}})
  tmp =  4*π*RADIUS^2/num_cells(model)
  sqrt(tmp)
end

function dx(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  horizontal = dx_horizontal(model)
  vertical = dx_vertical(model)
  horizontal*vertical
end

function dx_horizontal(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  horizontal = 4*π*RADIUS^2/(nc_horizontal(model)*6)
  sqrt(horizontal) ## quads so have to sqrt
end

function dx_vertical(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  vertical = THICKNESS/_nc_vertical(model)
  vertical ### single layer, so no sqrt
end

function nc(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  println("3D nc")
  nc_horizontal(model) + _nc_vertical(model)
end



# return square here so vertical is 'like' horitzontal
function nc_vertical(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  n = _nc_vertical(model)
  return Int(n^2)
end

# the actual number of cells in vertical per panel
function _nc_vertical(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  ncells_per_panel = nc_horizontal(model)
  n = num_cells(model)/6
  _n =  n /ncells_per_panel
  return Int(_n)
end
