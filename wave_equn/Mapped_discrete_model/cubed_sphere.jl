using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Test
include("structs.jl")

R = 1 # sphere radius
a = R/sqrt(3)

domain = (-a, a, -a, a)
partition = (10,10)
panel = CartesianDiscreteModel(domain,partition)
grid = get_grid(panel)
writevtk(Triangulation(panel),datadir("CubedSphere")*"/panel_og",append=false)

# function map_panel_to_latlong(xy)
#   x,y = xy
#   λ = atan(x,a)
#   θ = atan(y,a*sec(λ))
#   print(λ,θ)
#   return Point(λ,θ)
# end

# panel1_latlon_model = MappedDiscreteModel(panel,map_panel_to_latlong)
# writevtk(Triangulation(panel1_latlon_model),datadir("CubedSphere")*"/panel1_latlon",append=false)


function map_panel_to_sphere(x)
  # x,y = xy
  λ = atan(x[1],a)
  θ = atan(x[2],a*sec(λ))
  print(λ,θ)
  return VectorValue(R*cos(λ)*cos(θ), R*sin(λ)*cos(θ), R*sin(θ))
end


panel1_mapp_grid = MyMappedGrid(grid,map_panel_to_sphere)
panel1_cubed = MyMappedDiscreteModel(panel,panel1_mapp_grid)
writevtk(Triangulation(panel1_cubed),datadir("CubedSphere")*"/panel1_cubed",append=false)





function map_panel_xy_2_xyz(xy,panel)
  a,b=xy
  if panel==1
    x=Point(1.0,a,b)
  elseif panel==2
    x=Point(-a,1.0,b)
  elseif panel==3
    x=Point(-1.0,b,a)
  elseif panel==4
    x=Point(-b,-1.0,a)
  elseif panel==5
    x=Point(-b,a,1.0)
  elseif panel==6
    x=Point(-a,b,-1.0)
  end
  x
end


# map(x) = VectorValue(x[1],4*x[1],10*x[2])
