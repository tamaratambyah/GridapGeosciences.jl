using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../../../src/initialise.jl")

function bump_mats(a)

  # panel 1: front
  A1 = [0.0 1.0 0.0
      0.0 0.0 1.0]
  b1 = a .* [1.0
          0.0
          0.0]

  # panel 2: top
  A2 = [0.0 1.0 0.0
        -1.0 0.0 0.0]
  b2 = a .* [0.0
            0.0
            1.0]

  # panel 3: right
  A3 = [-1.0 0.0 0.0
        0.0 0.0 1.0]
  b3 = a .* [0.0
            1.0
            0.0]

  # panel 4: back
  A4 = [0.0 0.0 1.0
        0.0 1.0 0.0]
  b4 = -1.0 * b1

  # panel 5: bottom
  A5 = [-1.0 0.0 0.0
        0.0 1.0 0.0]
  b5 = -1.0 .* b2

  # panel 6: left
  A6 = [0.0 0.0 1.0
        -1.0 0.0 0.0]
  b6 = -1.0 .* b3

  As = [A1,A2,A3,A4,A5,A6]
  Bs = map(x->transpose(x),As)
  bs = [b1,b2,b3,b4,b5,b6]

  Abumps = map(x->TensorValue(x),As)
  Bbumps = map(x->TensorValue(x),Bs)
  bbumps = map(x->VectorValue(x),bs)

  return Abumps,Bbumps,bbumps
end


a = π/4
Abumps,Bbumps,bbumps = bump_mats(a)
RADIUS = 1.0*sqrt(3.0)
Panel2Sphere = [Panel123Sphere(RADIUS),
                Panel123Sphere(RADIUS),
                Panel123Sphere(RADIUS),
                Panel456Sphere(RADIUS),
                Panel456Sphere(RADIUS),
                Panel456Sphere(RADIUS)]

# rotation about Y axis by 3π/2
Ry = [0.0 0.0 -1.0
        0.0 1.0 0.0
        1.0 0.0 0.0]
# rotation about Z axis by π/2
Rz = [0.0 -1.0 0.0
      1.0 0.0 0.0
      0.0 0.0 1.0]
A_11 = Matrix{Float64}(I,3,3)

Rot_mats = [A_11,Ry,Rz,A_11,Ry,Rz]
Rs = map(x->TensorValue(x),Rot_mats)

## coarse cube
cube_grid_3D,topo_3D,face_labels_3D,panel_ids = coarse_cube_surface_3D(a)
cube_model_3D = UnstructuredDiscreteModel(cube_grid_3D,topo_3D,face_labels_3D)
writevtk(Triangulation(cube_model_3D),dir*"/ambient",append=false)

# cube_model_3D = Gridap.Adaptivity.refine(cube_model_3D)

## make parametric grid
panel_ids = get_panel_ids(cube_model_3D)
grid = get_grid(cube_model_3D)
topo = get_grid_topology(cube_model_3D)
cmaps = get_cell_map(grid)
cube_cell_coords = get_cell_coordinates(cube_model_3D)
## ambient grid
k = lazy_map(p-> PanelRotationField(Rs[p]) ∘ Panel2Sphere[p] ∘ BumpField(Abumps[p],Bbumps[p],bbumps[p]), panel_ids)
_cmaps = lazy_map(∘,k,cmaps)

cell_coords = lazy_map(evaluate,_cmaps,get_cell_ref_coordinates(grid))

cell_node_ids = get_cell_node_ids(grid)
nodes = similar(cell_coords, VectorValue{3,Float64}, num_nodes(grid))

for i in eachindex(cell_coords)
  ids = cell_node_ids[i]
  nodes[ids] .= cell_coords[i]
end

_grid = UnstructuredGrid(nodes,get_cell_node_ids(grid),get_reffes(grid),get_cell_type(grid),
  OrientationStyle(grid),nothing,_cmaps)

_topo = UnstructuredGridTopology(nodes,get_faces(topo,2,0),get_cell_type(topo),get_polytopes(topo),Gridap.Geometry.NonOriented())

model = UnstructuredDiscreteModel(_grid,_topo,FaceLabeling(_topo))

mask = map(x->x==1 || x ==2 || x == 3,panel_ids)
Ω =  Triangulation(model,mask)
writevtk(Ω,dir*"/ambient",append=false)



####
##


### parametric grid
k = lazy_map(p->  BumpField(Abumps[p],Bbumps[p],bbumps[p]), panel_ids)
_cmaps = lazy_map(∘,k,cmaps)

##
cell_coords = lazy_map(evaluate,_cmaps,get_cell_ref_coordinates(grid))


nodes = get_panel_1_nodes_from_coords(grid,cell_coords,panel_ids)

_grid = UnstructuredGrid(nodes,get_cell_node_ids(grid),get_reffes(grid),get_cell_type(grid),
  OrientationStyle(grid),nothing,_cmaps)

nodes_topo = get_panel_1_nodes_from_coords(topo,cell_coords,panel_ids)
_topo = UnstructuredGridTopology(nodes_topo,get_faces(topo,2,0),get_cell_type(topo),get_polytopes(topo),Gridap.Geometry.NonOriented())

model = UnstructuredDiscreteModel(_grid,_topo,FaceLabeling(_topo))



######### FE
function Gridap.CellData.get_cell_points(trian::Triangulation)
  println("my cell points")
  cell_ref_coords = get_cell_ref_coordinates(trian)
  cmaps = get_cell_map(trian)
  cell_phys_coords = lazy_map(evaluate,cmaps,cell_ref_coords)
  CellPoint(cell_ref_coords,cell_phys_coords,trian,ReferenceDomain())
end

Ω = Triangulation(model,mask)
dΩ = Measure(Ω,10)
pts = get_cell_points(Ω)

V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,3); conformity=:H1)
U = TrialFESpace(V)

uX_scalar(x) = x[1]*x[2]*x[3]
function u_scalar_ambient2parametric2(p::Int,uX::Function)
  function _u(xy)
    XYZ = Panel2Sphere[p](xy)
    uX(XYZ)
  end
end

cell_field = map(p->GenericField(u_scalar_ambient2parametric2(p,uX_scalar)),panel_ids)
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

for p in collect(1:6)
  mask = (panel_ids .== p)
  Ωp = Triangulation(model,mask)
  writevtk(Ωp,dir*"/u$p",cellfields=["u"=>ucf],append=false)
end

# _Jtcf = map(p->GenericField(∇(Panel2Sphere[p])),panel_ids)
# Jt = CellData.GenericCellField(_Jtcf,Ω,PhysicalDomain())
# J = Operation(transpose)(Jt)
# metric = Jt⋅J
# Jt(pts)[1]
# J(pts)[1]
# J(pts)[6]
# metric(pts)[1]
# metric(pts)[6]


function Jt1(x)
  TensorValue{2,3}(     RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*(sec(x[1]))^2,
                        RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*(sec(x[2]))^2,
                        RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                        RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*tan(x[2])*(sec(x[2]))^2,
                        RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*tan(x[1])*(sec(x[1]))^2,
                        RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *   (sec(x[1]))^2*(sec(x[2]))^2  )
end

function J1(x)
  TensorValue{3,2}(     RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*(sec(x[1]))^2,
                        RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                        RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*tan(x[1])*(sec(x[1]))^2,
                        RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*(sec(x[2]))^2,
                        RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*tan(x[2])*(sec(x[2]))^2,
                        RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *   (sec(x[1]))^2*(sec(x[2]))^2
                        )
end



function Jt2(x)
  TensorValue{2,3}(       RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*(sec(x[1]))^2,
                          RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*(sec(x[2]))^2,
                          RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*tan(x[2])*(sec(x[1]))^2,
                          RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                          RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                          RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*tan(x[1])*(sec(x[2]))^2  )
end

function J2(x)
  TensorValue{3,2}(       RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*(sec(x[1]))^2,
                          RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*tan(x[2])*(sec(x[1]))^2,
                          RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                          RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*(sec(x[2]))^2,
                          RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                          RADIUS/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*tan(x[1])*(sec(x[2]))^2
  )
end


function Jts(p)
  function _Jt(x)
    if p == 1 || p == 2 || p == 3
      return Jt1(x)
    else
      Jt2(x)
    end
  end
end

function Js(p)
  function _J(x)
    if p == 1 || p == 2 || p == 3
      return J1(x)
    else
      J2(x)
    end
  end
end

_metric1(x) = Jt1(x)⋅J1(x)
_metric2(x) = Jt2(x)⋅ J2(x)

function _invmetric1(p)
  function invm(x)
    if p == 1 || p == 2 || p == 3
      return   inv(_metric1(x))
    else
      return inv(_metric2(x))
    end
  end
end

function _sqrtmet(p)
  function sqm(x)
    if p == 1 || p == 2 || p == 2
    return sqrt(meas(_metric1(x)))
    else
      return sqrt(meas(_metric2(x)))
    end
  end
end

function panel_wise_cf(func::Function,Ω,panel_ids)
  _cf = map(p->GenericField(func(p)),panel_ids)
  return CellData.GenericCellField(_cf,Ω,PhysicalDomain())
end

Jt = panel_wise_cf(Jts,Ω,panel_ids)
J = panel_wise_cf(Js,Ω,panel_ids)

invM =  panel_wise_cf(_invmetric1,Ω,panel_ids)
sqrtM = panel_wise_cf(_sqrtmet,Ω,panel_ids)

pt = Point(π/4,π/4)


surflap =  1/sqrtM * divergence( sqrtM*( invM ⋅ gradient(ucf) ) )
surflap(pt)PanelPanel2Sphere2Sphere


rhs = ucf + surflap
sum(∫(rhs)dΩ)

mass(u,v) = ∫( (u*v)*sqrtM )dΩ
stiffnes(u,v) = ∫( ( (J ⋅invM ⋅ ∇(v)) ⋅ (J⋅ invM⋅ ∇(u) )) *sqrtM)dΩ
force(v) = ∫(  (rhs*v)*sqrtM )dΩ

helmholtz_biform(u,v) = mass(u,v) - stiffnes(u,v)
helmholtz_liform(v) = force(v)

op = AffineFEOperator(helmholtz_biform,helmholtz_liform,U,V)
uh = solve(LUSolver(),op)

e = l2(uh-ucf,dΩ)

for p in collect(1:6)
  mask = (panel_ids .== p)
  Ωp = Triangulation(model,mask)
  writevtk(Ωp,dir*"/u$p",cellfields=["u"=>ucf,"uh"=>uh,"rhs"=>rhs],append=false)
end
