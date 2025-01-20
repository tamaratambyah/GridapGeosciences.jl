using DrWatson
using Gridap
using Gridap.Geometry

function transfer_FE_map(cube_modelh::AdaptedDiscreteModel,mapH::FEFunction)

  map_basis = get_fe_basis(get_fe_space(mapH))
  map_trian = get_triangulation(map_basis)
  map_reffes = get_reffes(map_trian)
  order = get_order(map_reffes[1])

  T_vec = eltype(get_node_coordinates(cube_modelh))
  Vh = FESpace(cube_modelh,
              ReferenceFE(lagrangian,T_vec,order),
              conformity=:H1)

  maph = interpolate(mapH,Vh)

  return maph, order, Vh

end

# mapp(x) = Point(2*x[1],4*x[2],-1*x[3])
# cube_modelH = UnstructuredDiscreteModel(cube_grid)

# VH = FESpace(cube_modelH,
#             ReferenceFE(lagrangian,VectorValue{3,Float64},1),
#             conformity=:H1)
# mapH = interpolate(mapp,VH)
# cube_modelh = Gridap.Adaptivity.refine(cube_modelH)

# maph, order, Vh = transfer_FE_map(cube_modelh,mapH)
# _maph = interpolate(mapp,Vh)

# l2(w,dΩ) = sum( ∫( w⊙w )dΩ  )
# Ωh = Triangulation(cube_modelh)
# dΩh = Measure(Ωh,2*1+1)
# l2(maph - _maph,dΩh)
