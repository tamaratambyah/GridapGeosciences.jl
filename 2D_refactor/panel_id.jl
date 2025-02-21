using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Test
using LinearAlgebra
using FillArrays
using Plots
include("cube_surface_1_cell_per_panel.jl")


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()

cubemodel = UnstructuredDiscreteModel(cube_grid,topo,face_labels)
_cubemodelh = Gridap.Adaptivity.refine(cubemodel)
cubemodelh = Gridap.Adaptivity.refine(_cubemodelh)
cubemodelh2 = Gridap.Adaptivity.refine(cubemodelh)



function pid_ccam_2_pid_reg(pid_ccam::Int)
  if !(1 <= pid <= 6 )
    print("Invalid index CCAM panel indices are between 1 and 6")
    return
  end

  if pid_ccam == 1
    return pid_reg = 1
  elseif pid_ccam == 2
    return pid_reg = 5
  elseif pid_ccam == 3
    return pid_reg = 2
  elseif pid_ccam == 4
    return pid_reg = 3
  elseif pid_ccam == 5
    return pid_reg = 6
  elseif pid_ccam == 6
    return pid_reg = 4
  end


end

function get_xy_from_XYZ(pid_reg,X,Y,Z)
  # convert from cube to local panel coordinates using Ronchi1996 equns (6)-(13)

  # Note pid_reg is the panel index in Ronchi1996
  # Note X,Y,Z are cartesian coordinates of cube (Ronchi1996 uses x,y,z)
  # Note x,y are local to the panel (Ronchi1996 uses X,Y)


  if pid_reg == 1
    x = Y./X
    y = Z./X
  elseif pid_reg == 2
    x = -X./Y
    y = Z./Y
  elseif pid_reg == 3
    x = Y./X
    y = -Z./X
  elseif pid_reg == 4
    x = -X./Y
    y = -Z./Y
  elseif pid_reg == 5
    x = Y./Z
    y = - X./Z
  elseif pid_reg == 6
    x = -Y./Z
    y = -X./Z
  end

  return x,y
end

function get_panel_info(X,Y,Z)
  panel_assigned = true

  if length(findall(i->i==1,X)) == length(X)
    println("ccam panel1")

    ccam_pidx = 1
    og_pidx = 1

    λϕc = Point(0.0,0.0)

  elseif length(findall(i->i==-1,X)) == length(X)
    println("ccam panel4")

    ccam_pidx = 4
    og_pidx = 3

    λϕc = Point( (og_pidx-1)*π/2,0.0)

  elseif length(findall(i->i==1,Y)) == length(Y)
    println("ccam panel3")
    ccam_pidx = 3
    og_pidx = 2

    λϕc = Point( (og_pidx-1)*π/2,0.0)

  elseif length(findall(i->i==-1,Y)) == length(Y)
    println("ccam panel6")

    ccam_pidx = 6
    og_pidx = 4

    λϕc = Point( (og_pidx-1)*π/2,0.0)

  elseif length(findall(i->i==1,Z)) == length(Z)
    println("ccam panel2")
    ccam_pidx = 2
    og_pidx = 5

    λϕc = Point(0.0,π/2)

  elseif length(findall(i->i==-1,Z)) == length(Z)
    println("ccam panel5")
    ccam_pidx = 5
    og_pidx = 6

    λϕc = Point(0.0,-π/2)

  else
    println("not assigned to a panel")
    panel_assigned = false
  end

  return panel_assigned, ccam_pidx, og_pidx, λϕc

end


function get_XYZ_from_cell(cell)
  X = map(x->x[1],cell)
  Y = map(x->x[2],cell)
  Z = map(x->x[3],cell)
  return X,Y,Z
end


function compute_latlon(x,y,λϕc)

  # convert to angular variables ζ,η (using Rochni1996 eq (1)-(2) )
  ζ = atan.(x)
  η = atan.(y)

  # convert to lat-lon using Girdaldo2003 eq (32)
  # Note, using notation from Girdalod2003 where λ = lon, ϕ = lat
  λg = ζ
  ϕg = asin.( y./ ( (ones(size(x))+ (x).^2 + (y).^2 ).^(0.5) ) )

  λc,ϕc = λϕc # get centroid of faces from panel id

  # rotate using Girdaldo2003 eq (33)
  ϕ = asin.( sin.(ϕg)*cos(ϕc) + cos.(ϕg).*cos.(λg)*sin(ϕc)      )
  _λ = λc*ones(size(λg)) + atan.( cos.(ϕ).*sin.(λg), cos.(ϕg).*cos.(λg)*cos(ϕc) - sin.(ϕg)*sin(ϕc))
  λ = rem2pi.(_λ,RoundNearest) # put in interval [-π,π]

  return λ,ϕ
end

function map_cube_cell_2_latlon(cell)
  X,Y,Z = get_XYZ_from_cell(cell)
  panel_assigned, ccam_pidx, og_pidx, λϕc = get_panel_info(X,Y,Z)
  @assert panel_assigned "Panel not assigned!"
  x,y = get_xy_from_XYZ(og_pidx,X,Y,Z)
  λ,ϕ = compute_latlon(x,y,λϕc)
  return λ,ϕ,ccam_pidx
end

### Check for node4 = (1,1,1) -> cells 1,2,3
model = cubemodel
cell_phys_coords = get_cell_coordinates(get_grid(model))./1

cell = cell_phys_coords[1]
λ,ϕ,ccam_pidx = map_cube_cell_2_latlon(cell)
λ1,ϕ1 = λ[4],ϕ[4]

cell = cell_phys_coords[2]
λ,ϕ,ccam_pidx = map_cube_cell_2_latlon(cell)
λ2,ϕ2 = λ[2],ϕ[2]

cell = cell_phys_coords[3]
λ,ϕ,ccam_pidx = map_cube_cell_2_latlon(cell)
λ3,ϕ3 = λ[1],ϕ[1]

λ1-λ2
λ1-λ3
λ2-λ3

ϕ1-ϕ2
ϕ1-ϕ3
ϕ2-ϕ3


### Check for node7 = (-1,1,-1) -> cells 3,4,5
cell = cell_phys_coords[3]
λ,ϕ,ccam_pidx = map_cube_cell_2_latlon(cell)
λ1,ϕ1 = λ[4],ϕ[4]

cell = cell_phys_coords[4]
λ,ϕ,ccam_pidx = map_cube_cell_2_latlon(cell)
λ2,ϕ2 = λ[2],ϕ[2]

cell = cell_phys_coords[5]
λ,ϕ,ccam_pidx = map_cube_cell_2_latlon(cell)
λ3,ϕ3 = λ[1],ϕ[1]

λ1-λ2
λ1-λ3
λ2-λ3

ϕ1-ϕ2
ϕ1-ϕ3
ϕ2-ϕ3


### Check for node1 = (1,-1,-1) -> cells 1,5,6
cell = cell_phys_coords[1]
λ,ϕ,ccam_pidx = map_cube_cell_2_latlon(cell)
λ1,ϕ1 = λ[1],ϕ[1]

cell = cell_phys_coords[5]
λ,ϕ,ccam_pidx = map_cube_cell_2_latlon(cell)
λ2,ϕ2 = λ[4],ϕ[4]

cell = cell_phys_coords[6]
λ,ϕ,ccam_pidx = map_cube_cell_2_latlon(cell)
λ3,ϕ3 = λ[2],ϕ[2]

λ1-λ2
λ1-λ3
λ2-λ3

ϕ1-ϕ2
ϕ1-ϕ3
ϕ2-ϕ3




### map all cells
model = cubemodelh
cell_phys_coords = get_cell_coordinates(get_grid(model))./1
Ncells = length(cell_phys_coords)

λs = [Vector{Float64}(undef,4)]
ϕs = [Vector{Float64}(undef,4)]

Xs = [Vector{Float64}(undef,4)]
Ys = [Vector{Float64}(undef,4)]
Zs = [Vector{Float64}(undef,4)]

Xc = [Vector{Float64}(undef,4)]
Yc = [Vector{Float64}(undef,4)]
Zc = [Vector{Float64}(undef,4)]

pids = Vector{Int32}()

for i = 1:Ncells

  cell = cell_phys_coords[i]
  X,Y,Z = get_XYZ_from_cell(cell)
  λ,ϕ,ccam_pidx = map_cube_cell_2_latlon(cell)

  # get 3D sphere points
  Xsphere = cos.(λ).*cos.(ϕ)
  Ysphere = sin.(λ).*cos.(ϕ)
  Zsphere = sin.(ϕ)

  push!(λs,λ)
  push!(ϕs,ϕ)
  push!(Xs,Xsphere)
  push!(Ys,Ysphere)
  push!(Zs,Zsphere)
  push!(Xc,X)
  push!(Yc,Y)
  push!(Zc,Z)
  push!(pids,ccam_pidx)
end


### plotting
markers = [:circle, :rect, :diamond,  :utriangle,  :x, :cross,  ]
colors = palette(:tab10)



plot()
for i in 1:length(pids)
  p = pids[i]
  scatter!(λs[i+1,:],ϕs[i+1,:],markershape=markers[p],mc=colors[p])
  # i+1 index to remove first row
end
plot!(show=true,legend=false)
savefig("lat-lon.png")

plot()
for i in 1:length(pids)
  p = pids[i]
  scatter!(Xs[i+1,:],Ys[i+1,:],Zs[i+1,:],
      markershape=markers[p],mc=colors[p])
  # scatter!(Xc[i+1,:],Yc[i+1,:],Zc[i+1,:],
  #     markershape=markers[p],mc=colors[p])
  # i+1 index to remove first row
end
plot!(show=true,legend=false)
savefig("sphere.png")


# for i = 1:Ncells
# # i = 1
#   panel_assigned = true
#   cell = cell_phys_coords[i]
#   X = map(x->x[1],cell)
#   Y = map(x->x[2],cell)
#   Z = map(x->x[3],cell)

#   if length(findall(i->i==1,X)) == length(X)
#     println("ccam panel1")
#     # a = Y
#     # b = Z

#     ccam_pidx = 1
#     og_pidx = 1

#     λϕc = Point( (og_pidx-1)*π/2,0.0)

#   elseif length(findall(i->i==-1,X)) == length(X)
#     println("ccam panel4")
#     # a = -Z
#     # b = -Y

#     ccam_pidx = 4
#     og_pidx = 3

#     λϕc = Point( (og_pidx-1)*π/2,0.0)

#   elseif length(findall(i->i==1,Y)) == length(Y)
#     println("ccam panel3")
#     # a = -Z
#     # b = -X

#     ccam_pidx = 3
#     og_pidx = 2

#     λϕc = Point( (og_pidx-1)*π/2,0.0)

#   elseif length(findall(i->i==-1,Y)) == length(Y)
#     println("ccam panel6")
#     # a = X
#     # b = Z

#     ccam_pidx = 6
#     og_pidx = 4

#     λϕc = Point( (og_pidx-1)*π/2,0.0)

#   elseif length(findall(i->i==1,Z)) == length(Z)
#     println("ccam panel2")
#     # a = Y
#     # b = -X

#     ccam_pidx = 2
#     og_pidx = 5

#     λϕc = Point(0.0,π/2)

#   elseif length(findall(i->i==-1,Z)) == length(Z)
#     println("ccam panel5")
#     # a = X
#     # b = -Y

#     ccam_pidx = 5
#     og_pidx = 6

#     λϕc = Point(0.0,-π/2)

#   else
#     println("not assigned to a panel")
#     panel_assigned = false
#   end

#   # λc,ϕc = λϕc
#   # ϕg = asin.( b./ ( (ones(size(a))+ (a).^2 + (b).^2 ).^(0.5) ) )

#   # sinϕg = b./ ( (ones(size(a))+ (a).^2 + (b).^2 ).^(0.5) )
#   # cos.(ϕg).*cos.(λg)*sin(ϕc)


#   @assert panel_assigned "Panel not assigned!"

#   α = atan.(a)
#   β = atan.(b)

#   λg = α
#   ϕg = asin.( b./ ( (ones(size(a))+ (a).^2 + (b).^2 ).^(0.5) ) )

#   λc,ϕc = λϕc
#   ϕ = asin.( sin.(ϕg)*cos(ϕc) + cos.(ϕg).*cos.(λg)*sin(ϕc)      )
#   λ = λc*ones(size(λg)) + atan.( cos.(ϕ).*sin.(λg), cos.(ϕg).*cos.(λg)*cos(ϕc) - sin.(ϕg)*sin(ϕc))

#   # λ = rem2pi.(_λ,RoundNearest) # put in interval [-π,π]


#   Xsphere = cos.(λ).*cos.(ϕ)
#   Ysphere = sin.(λ).*cos.(ϕ)
#   Zsphere = sin.(ϕ)

#   push!(λs,λ)
#   push!(ϕs,ϕ)
#   push!(Xs,Xsphere)
#   push!(Ys,Ysphere)
#   push!(Zs,Zsphere)
#   push!(Xc,X)
#   push!(Yc,Y)
#   push!(Zc,Z)
#   push!(pids,ccam_pidx)
# end

lon =  collect(Iterators.flatten(λs[2:end]))
lat = collect(Iterators.flatten(ϕs[2:end]))

all = [lon lat]

ix = unique(i -> all[i, :], axes(all, 1))
ll_pairs = unique(all,dims=1)

l = unique(x-> round(x,digits=16),lon)
plot()
scatter!(1:length(l),l)

k = unique(x-> round(x,digits=17),lat)
plot()
scatter!(1:length(k),k)
