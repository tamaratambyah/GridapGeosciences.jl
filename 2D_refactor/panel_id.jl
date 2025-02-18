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
include("cube_surface_1_cell_per_panel.jl")


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()

cubemodel = UnstructuredDiscreteModel(cube_grid,topo,face_labels)
# _cubemodelh = Gridap.Adaptivity.refine(cubemodel)
# cubemodelh = Gridap.Adaptivity.refine(_cubemodelh)

model = cubemodel
cell_phys_coords = get_cell_coordinates(get_grid(model))./1

using Plots
plot()
# for i = 1:length(cell_phys_coords)
i = 1
  cell = cell_phys_coords[i]
  X = map(x->x[1],cell)
  Y = map(x->x[2],cell)
  Z = map(x->x[3],cell)

  if length(findall(i->i==1,X)) == length(X)
    println("panel1")
    a = Y
    b = Z
    λϕc = Point(0,0)
  elseif length(findall(i->i==-1,X)) == length(X)
    println("panel4")
    a = -Z
    b = -Y
    λϕc = Point(π/2,0)
  elseif length(findall(i->i==1,Y)) == length(Y)
    println("panel3")
    a = -Z
    b = -X
    λϕc = Point(π,0)
  elseif length(findall(i->i==-1,Y)) == length(Y)
    println("panel6")
    a = X
    b = Z
    λϕc = Point(3*π/2,0)
  elseif length(findall(i->i==1,Z)) == length(Z)
    println("panel2")
    a = Y
    b = -X
    λϕc = Point(0,π/2)
  elseif length(findall(i->i==-1,Z)) == length(Z)
    println("panel5")
    a = X
    b = -Y
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
  λ = λc*ones(size(λg)) + atan.( cos.(ϕ).*sin.(λg), cos.(ϕg).*cos.(λg)*cos(ϕc) - sin.(ϕg)*sin(ϕc)     )


  scatter!(λ,ϕ)

# end
plot!(show=true)
