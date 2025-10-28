using Gridap
using GridapDistributed
using GridapGeosciences
using Gridap.Fields
using Gridap.Helpers

# radius(XYZ) = sqrt(XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2] + XYZ[3]*XYZ[3])

# RADIUS = 1.0 # radius of inner sphere (hardcoded at moment)
# RADIUS_OUTER = 2.0 # radius of outer sphere

function forward_map_3D(p::Int,γαβ)
  @check length(γαβ) == 3 "\n Not 3D point"

  #### recall the first coordinate in P6est is the extrusion!
  γ,α,β = γαβ
  # α,β,γ = γαβ

  #### compute XYZ point on surface of inner sphere using 2D forward_map
  αβ = Point(α,β)
  XYZ_surf = forward_map_2D(p,αβ)

  # radius_surf = radius(XYZ)
  # radius_surf = RADIUS

  #### extrude surface point in radial direction
  # XYZ_surf + (RADIUS_OUTER-radius_surf)*γ*normal_vec(XYZ_surf)
  return XYZ_surf + γ* normal_vec(XYZ_surf)


end

forward_map_3D(p::Int) = γαβ -> forward_map_3D(p,γαβ)
forward_jacobian_3D(p::Int,γαβ) = transpose( gradient(forward_map_3D(p))(γαβ) )
forward_jacobian_3D(p::Int) = γαβ -> forward_jacobian_3D(p,γαβ)

# function forward_jacobian_3D(p::Int,γαβ)

#   γ,α,β = γαβ
#   rho = sqrt(1 + (tan(α))^2 + (tan(β))^2 )
#   drho_da = - tan(α)*(sec(α))^2 / ( rho^3 )
#   drho_db = - tan(β)*(sec(β))^2 / ( rho^3 )

#   if p == 1
#     # X = 1/rho
#     # Y = 1/rho * tan(α)
#     # Z = 1/rho * tan(β)

#     dXda = drho_da
#     dXdb = drho_db
#     dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
#     dYdb = drho_db*tan(α)
#     dZda = drho_da*tan(β)
#     dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2

#     dXdg = 1/RADIUS*( 1/rho )
#     dYdg = 1/RADIUS*( 1/rho * tan(α) )
#     dZdg = 1/RADIUS*( 1/rho * tan(β) )
#   elseif p == 2
#     # X = -1/rho * tan(β)
#     # Y = 1/rho * tan(α)
#     # Z = 1/rho

#     dXda = -( drho_da*tan(β) )
#     dXdb = -( drho_db*tan(β) + 1/rho*(sec(β))^2 )
#     dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
#     dYdb = drho_db*tan(α)
#     dZda = drho_da
#     dZdb = drho_db

#     dXdg = 1/RADIUS*( -1/rho * tan(β) )
#     dYdg = 1/RADIUS*( 1/rho * tan(α) )
#     dZdg = 1/RADIUS*( 1/rho )
#   elseif p == 3
#     # X = -1/rho * tan(α)
#     # Y = 1/rho
#     # Z = 1/rho * tan(β)

#     dXda = -( drho_da*tan(α) + 1/rho*(sec(α))^2  )
#     dXdb = -( drho_db*tan(α) )
#     dYda = drho_da
#     dYdb = drho_db
#     dZda = drho_da*tan(β)
#     dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2

#     dXdg = 1/RADIUS*( -1/rho * tan(α) )
#     dYdg = 1/RADIUS*( 1/rho )
#     dZdg = 1/RADIUS*( 1/rho * tan(β) )
#   elseif p == 4
#     # X = -1/rho
#     # Y = 1/rho * tan(β)
#     # Z = 1/rho * tan(α)

#     dXda = -drho_da
#     dXdb = -drho_db
#     dYda = drho_da*tan(β)
#     dYdb = drho_db*tan(β) + 1/rho*(sec(β))^2
#     dZda = drho_da*tan(α) + 1/rho*(sec(α))^2
#     dZdb = drho_db*tan(α)

#     dXdg = 1/RADIUS*( -1/rho )
#     dYdg = 1/RADIUS*( 1/rho * tan(β) )
#     dZdg = 1/RADIUS*( 1/rho * tan(α) )
#   elseif p == 5
#     # X = -1/rho * tan(α)
#     # Y = 1/rho * tan(β)
#     # Z = -1/rho

#     dXda = -(drho_da*tan(α) + 1/rho*(sec(α))^2)
#     dXdb = -(drho_db*tan(α) )
#     dYda = drho_da*tan(β)
#     dYdb = drho_db*tan(β) + 1/rho*(sec(β))^2
#     dZda = -drho_da
#     dZdb = -drho_db

#     dXdg = 1/RADIUS*( -1/rho * tan(α) )
#     dYdg = 1/RADIUS*( 1/rho * tan(β) )
#     dZdg = 1/RADIUS*( -1/rho )
#   elseif p == 6
#     # X = -1/rho * tan(β)
#     # Y = -1/rho
#     # Z = 1/rho * tan(α)

#     dXda = -(drho_da*tan(β))
#     dXdb = -(drho_db*tan(β) + 1/rho*(sec(β))^2)
#     dYda = -drho_da
#     dYdb = -drho_db
#     dZda = drho_da*tan(α) + 1/rho*(sec(α))^2
#     dZdb = drho_db*tan(α)

#     dXdg = 1/RADIUS*( -1/rho * tan(β) )
#     dYdg = 1/RADIUS*( -1/rho )
#     dZdg = 1/RADIUS*( 1/rho * tan(α) )
#   end

#   ## J = [dXdg dXda dXdb
#   ##      dYdg dYda dYdb
#   ##      dZdg dZda dZdb  ]
#   ## As a TensorValue data = (dXdg,dYdg,dZdg, dXda,dYda,dZda, dXdb,dYdb,dZdb)

#   RADIUS*TensorValue{3,3}(dXdg,dYdg,dZdg, (1+γ)*dXda, (1+γ)*dYda, (1+γ)*dZda, (1+γ)*dXdb,(1+γ)*dYdb,(1+γ)*dZdb)

# end
