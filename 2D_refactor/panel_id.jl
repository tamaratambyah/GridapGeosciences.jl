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
# using PlotlyJS
using Plots
include("cube_surface_1_cell_per_panel.jl")


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()

cubemodel = UnstructuredDiscreteModel(cube_grid,topo,face_labels)
_cubemodelh = Gridap.Adaptivity.refine(cubemodel)
cubemodelh = Gridap.Adaptivity.refine(_cubemodelh)
cubemodelh2 = Gridap.Adaptivity.refine(cubemodelh)

model = cubemodelh
cell_phys_coords = get_cell_coordinates(get_grid(model))./1

markers = [:circle, :rect, :diamond,  :utriangle,  :x, :cross,  ]
colors = palette(:tab10)

plot()
for i = 1:length(cell_phys_coords)
  cell = cell_phys_coords[i]
  X = map(x->x[1],cell)
  Y = map(x->x[2],cell)
  Z = map(x->x[3],cell)

  if length(findall(i->i==1,X)) == length(X)
    println("ccam panel1")
    a = Y
    b = Z

    _a = a
    _b = b

    ccam_pidx = 1
    og_pidx = 1

    λϕc = Point( (og_pidx-1)*π/2,0)

  elseif length(findall(i->i==-1,X)) == length(X)
    println("ccam panel4")
    a = -Z
    b = -Y

    _a = -a
    _b = b

    ccam_pidx = 4
    og_pidx = 3

    λϕc = Point( (og_pidx-1)*π/2,0)

  elseif length(findall(i->i==1,Y)) == length(Y)
    println("ccam panel3")
    a = -Z
    b = -X

    _a = a
    _b = b

    ccam_pidx = 3
    og_pidx = 2

    λϕc = Point( (og_pidx-1)*π/2,0)

  elseif length(findall(i->i==-1,Y)) == length(Y)
    println("ccam panel6")
    a = X
    b = Z

    _a = a
    _b = b

    ccam_pidx = 6
    og_pidx = 4

    λϕc = Point( (og_pidx-1)*π/2,0)

  elseif length(findall(i->i==1,Z)) == length(Z)
    println("ccam panel2")
    a = Y
    b = -X

    _a = a
    _b = b

    ccam_pidx = 2
    og_pidx = 5

    λϕc = Point(0,π/2)

  elseif length(findall(i->i==-1,Z)) == length(Z)
    println("ccam panel5")
    a = X
    b = -Y

    _a = -a
    _b = b

    ccam_pidx = 5
    og_pidx = 6

    λϕc = Point(0,-π/2)

  else
    println("not assigned to a panel")
  end

  α = atan.(a)
  β = atan.(b)

  λg = α
  ϕg = asin.( b./ ( (ones(size(a))+ (a).^2 + (b).^2 ).^(0.5) ) )

  λc,ϕc = λϕc
  ϕ = asin.( sin.(ϕg)*cos(ϕc) + cos.(ϕg).*cos.(λg)*sin(ϕc)      )
  _λ = λc*ones(size(λg)) + atan.( cos.(ϕ).*sin.(λg), cos.(ϕg).*cos.(λg)*cos(ϕc) - sin.(ϕg)*sin(ϕc)     )

  λ = rem2pi.(_λ,RoundNearest) # put in interval [-π,π]

  # ϕ = ϕg + ϕc*ones(size(ϕg))
  # λ = λg + λc*ones(size(λg))

  Xsphere = cos.(λ).*cos.(ϕ)
  Ysphere = sin.(λ).*cos.(ϕ)
  Zsphere = sin.(ϕ)
  scatter!(λ,ϕ,markershape=markers[ccam_pidx],mc=colors[ccam_pidx])
  # scatter!(Xsphere,Ysphere,Zsphere;c=colors[ccam_pidx],s=1)
end
plot!(show=true,legend=false)
savefig("sphere.png")
