using DrWatson
using Gridap
using GridapDistributed
using GridapP4est
using GridapGeosciences
using GridapGeosciences.Distributed
using Gridap.Helpers
using MPI
using PartitionedArrays

using DrWatson
dir = datadir("Mountain")
!isdir(dir) && mkdir(dir)

function height_func(θ,ϕ)
  ζ = 0.0

  a_e = 6.37e6 # m
  g = 9.8 # m/2
  ω = 7.29e-5 #s^-1
  H_0 = 5960.0 #2.94e4/g #m
  u_0 = 20 #2*π*a_e/T #m/s

  L = a_e
  _τ = 1/ω#

  _g = g*_τ^2/L
  _ω = ω*_τ
  _H_0 = H_0/L
  _u0 = u_0/L*_τ

  h  = -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ)
  _H_0 - (_ω*_u0 + 0.5*_u0*_u0)*h*h/_g

end

function height_func(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  height_func(θ,ϕ)
end

# Topography
function topography_func(θ,ϕ)

  b_0 = 2000.0 #m
  L = 6.37e6 # m
  _b0 = b_0/L

  θc    =  -π/2 #0.0
  ϕc    =  π/6
  R   = π/9.0
  rsq   = (ϕ - ϕc)*(ϕ - ϕc) + (θ - θc)*(θ - θc)
  r2 = min(R^2,rsq)

  r = sqrt(r2)
  b = _b0*(1.0 - r/R)
  b
end

function topography_func(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  topography_func(rem2pi(θ,RoundNearest),rem2pi(ϕ,RoundNearest))
end

# mountain(θ,ϕ) =  topography_func(θ,ϕ)# height_func(θ,ϕ) +

θs = LinRange(0, 2*π,  40)
ϕs = LinRange(-π/2, π/2, 40)


# using GLMakie


### plot topography_func
zs = [topography_func(θ,ϕ) for θ in θs, ϕ in ϕs]

f = Figure()
ax = Axis3(f[1, 1], zlabel = "height", xlabel = "lon", ylabel="lat")
s = surface!(ax,θs,ϕs,zs)
Colorbar(f[1, 2], s)
save(dir*"/topograph.png", f)

f = Figure()
ax = Axis(f[1, 1],xlabel="lon",ylabel="lat")
co = contourf!(ax,θs,ϕs,zs,levels=10,
    extendlow = :auto, extendhigh = :auto)
Colorbar(f[1, 2], co)
save(dir*"/topograph_contourf.png", f)


f = Figure()
ax = Axis(f[1, 1],xlabel="lon",ylabel="lat")
co = contour!(ax,θs,ϕs,zs,levels=10)
save(dir*"/topograph_contour.png", f)

### plot initial free surface
zs = [height_func(θ,ϕ) for θ in θs, ϕ in ϕs]

f = Figure()
ax = Axis3(f[1, 1], zlabel = "free surface", xlabel = "lon", ylabel="lat")
s = surface!(ax,θs,ϕs,zs)
Colorbar(f[1, 2], s)
save(dir*"/free_surface.png", f)

f = Figure()
ax = Axis(f[1, 1],xlabel="lon",ylabel="lat")
co = contourf!(ax,θs,ϕs,zs,levels=10,
    extendlow = :auto, extendhigh = :auto)
Colorbar(f[1, 2], co)
save(dir*"/free_surface_contourf.png", f)


f = Figure()
ax = Axis(f[1, 1],xlabel="lon",ylabel="lat")
co = contour!(ax,θs,ϕs,zs,levels=10)
save(dir*"/free_surface_contour.png", f)

### plot initial depth
zs = [(height_func(θ,ϕ)-topography_func(θ,ϕ)) for θ in θs, ϕ in ϕs]

f = Figure()
ax = Axis3(f[1, 1], zlabel = "depth", xlabel = "lon", ylabel="lat")
s = surface!(ax,θs,ϕs,zs)
Colorbar(f[1, 2], s)
save(dir*"/depth.png", f)

f = Figure()
ax = Axis(f[1, 1],xlabel="lon",ylabel="lat")
co = contourf!(ax,θs,ϕs,zs,levels=10,
    extendlow = :auto, extendhigh = :auto)
Colorbar(f[1, 2], co)
save(dir*"/depth_contourf.png", f)


f = Figure()
ax = Axis(f[1, 1],xlabel="lon",ylabel="lat")
co = contour!(ax,θs,ϕs,zs,levels=10)
save(dir*"/depth_contour.png", f)


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

num_horizontal_uniform_refinements = 5
num_vertical_uniform_refinements = 3
omodel = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
	num_horizontal_uniform_refinements=num_horizontal_uniform_refinements,
  num_vertical_uniform_refinements=num_vertical_uniform_refinements);




function GridapGeosciences.Helpers.forward_map_3D(p::Int,γαβ)


  println("mountain map")
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
  ht = panel_to_cartesian(topography_func)(p)(αβ)

  return XYZ_surf + (γ* normal_vec(XYZ_surf))*ht


end

dmodel = omodel.parametric_dmodel
Ω_panel = Triangulation(dmodel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/extruded_model",append=false,geo_map=cell_geo_map)
