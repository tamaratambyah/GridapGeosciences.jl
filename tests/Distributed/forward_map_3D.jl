using Gridap
using GridapDistributed
using GridapGeosciences
using Gridap.Fields
using Gridap.Helpers

radius(XYZ) = sqrt(XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2] + XYZ[3]*XYZ[3])

RADIUS = 1.0 # radius of inner sphere (hardcoded at moment)
RADIUS_OUTER = 2.0 # radius of outer sphere

function forward_map_3D(p::Int,γαβ)
  @check length(γαβ) == 3 "\n Not 3D point"

  #### recall the first coordinate in P6est is the extrusion!
  γ,α,β = γαβ

  #### compute XYZ point on surface of inner sphere using 2D forward_map
  αβ = Point(α,β)
  XYZ_surf = forward_map(p,αβ)

  # radius_surf = radius(XYZ)
  radius_surf = RADIUS

  #### extrude surface point in radial direction
  XYZ_surf #+ (RADIUS_OUTER-radius_surf)*γ*normal_vec(XYZ_surf)

end

A_3D_2_2D_panel = TensorValue{2,3}(0.0,0.0, 1.0,0.0, 0.0,1.0)

γαβ = Point(1.0,  π/4,π/4)
αβ = A_3D_2_2D_panel ⋅ γαβ
forward_map_3D(1,γαβ)
forward_map(1,αβ)


function geo_map_func_3D(trian::GridapDistributed.DistributedTriangulation)
  println("distributed 3D geo map")

  model = get_background_model(trian)
  # owned_panel_ids = get_owned_panel_ids(model)
  owned_panel_ids = get_panel_ids(trian)
  return geo_map_func_3D(owned_panel_ids)
end


function geo_map_func_3D(owned_panel_ids::AbstractArray)
  println("distributed 3D geo map")

  @assert typeof(owned_panel_ids) <: DebugArray || typeof(owned_panel_ids) <: MPIArray "\n Not distributed panel ids"

  cell_geo_map = map(owned_panel_ids) do pid
    return lazy_map(p -> ForwardMapPanel3D(p) , pid)
  end
  return cell_geo_map
end

function geo_map_func_3D(panel_ids::AbstractArray{Int})
  println("serial 3D geo map")
  return lazy_map(p -> ForwardMapPanel3D(p), panel_ids)
end


struct ForwardMapPanel3D{A}  <: Field
  p::A
end

function Gridap.Arrays.return_cache(f::ForwardMapPanel3D,cellx::AbstractArray{<:VectorValue{3}})
  y = similar(cellx,VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::ForwardMapPanel3D,cellx::AbstractArray{<:VectorValue{3}} )
  p = f.p
  y = cache
  map!(x-> forward_map_3D(p,x),
      y, cellx  )
  return y
end


function Gridap.Arrays.return_cache(f::ForwardMapPanel3D,x::VectorValue{3})
  y = zero(VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::ForwardMapPanel3D,x::VectorValue{3})
  p = f.p
  y = cache
  y = forward_map_3D(p,x)
  return y
end
