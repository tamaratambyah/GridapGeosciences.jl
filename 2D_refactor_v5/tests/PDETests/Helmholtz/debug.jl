using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../../../src/initialise.jl")



RADIUS = 1.0
a = RADIUS/sqrt(3)
cube_grid_3D,topo_3D,face_labels_3D,panel_ids = coarse_cube_surface_3D(a)
cube_model_3D = UnstructuredDiscreteModel(cube_grid_3D,topo_3D,face_labels_3D)
A_bump, B_bump, b_bump = bump_matrics(a)


### apply refinement
cube_model_3D = Gridap.Adaptivity.refine(cube_model_3D)
# get_cell_coordinates(cube_model_3D)
panel_ids = get_panel_ids(cube_model_3D)

#### make 2D surface
cube_grid_3D = get_grid(cube_model_3D)

cmaps_3D = get_cell_map(cube_grid_3D)

include("../../../src/Fields/Tan.jl")
k = lazy_map(p-> TanField(a) ∘ BumpField(A_bump,B_bump,b_bump) ∘ PanelRotationField(rp1_3D[p]), panel_ids)

parametric_cell_map = lazy_map(∘,k,cmaps_3D)

cell_coords_3D = get_cell_coordinates(cube_grid_3D)
cell_coords_panel1_2D = get_cube_nodes_2D(cell_coords_3D,panel_ids)
nodes_2D = get_panel_1_nodes_from_coords(cube_grid_3D,cell_coords_panel1_2D,panel_ids)


# nodes_2D = zeros(VectorValue{2,Float64}, num_nodes(cube_grid_3D))

topo_2D = UnstructuredGridTopology(nodes_2D,
    get_cell_node_ids(cube_grid_3D),get_cell_type(cube_grid_3D),[QUAD],Gridap.Geometry.NonOriented())
face_labels_2D = get_face_labeling(cube_model_3D)

cube_grid_2D = Gridap.Geometry.UnstructuredGrid(nodes_2D,
      get_cell_node_ids(cube_grid_3D),get_reffes(cube_grid_3D),get_cell_type(cube_grid_3D),Gridap.Geometry.NonOriented(),
      nothing,parametric_cell_map)

cube_model_2D = UnstructuredDiscreteModel(cube_grid_2D,topo_2D,face_labels_2D)
writevtk(cube_model_2D,dir*"/panel",append=false)

evaluate(parametric_cell_map[1],Point(0,0))

# cc = get_cell_map(Ω)
# evaluate(cc[1],Point(0,0))

### analytic
function uαβ_scalar(p)
  function _h(αβ)
    if  p == 2 || p == 3 || p == 4
      return -αβ[1]*αβ[2]
    else
      return αβ[1]*αβ[2]
    end
  end
end

function lapαβ_scalar(p)
  function _h(αβ)
    _sq_meas_func = sq_meas_func(cubedsphere)
    _inv_metric_func = inv_metric_func(cubedsphere)
    if p == 2 || p == 3 || p == 4
      grad(αβ) = VectorValue(-αβ[2],-αβ[1])
      sg(αβ) = _inv_metric_func(αβ)⋅grad(αβ)
      return 1/_sq_meas_func(αβ) * ( _sq_meas_func(αβ) ⋅ (divergence( sg)(αβ))
      + gradient( _sq_meas_func )(αβ)⋅ ( sg(αβ)  )
      )
    else
      _grad(αβ) = VectorValue(αβ[2],αβ[1])
      _sg(αβ) = _inv_metric_func(αβ)⋅_grad(αβ)
      return 1/_sq_meas_func(αβ) * ( _sq_meas_func(αβ) ⋅ (divergence( _sg)(αβ))
              + gradient( _sq_meas_func )(αβ)⋅ ( _sg(αβ)  )
      )
    end
  end
end



##### FE problem

Ω = Triangulation(cube_model_2D)
Γ = BoundaryTriangulation(cube_model_2D;tags="boundary")

writevtk(Γ,dir*"/panel",append=false)


m = Metric(cubedsphere,Ω)
mΓ = Metric(cubedsphere,Γ)

dΓ = Measure(Γ,10)
n_Γ = get_normal_vector(Γ)

dΩ = Measure(Ω,10)

pts = get_cell_points(Ω)

cell_field = map(p->GenericField(uαβ_scalar(p)),panel_ids)
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

lap_ucf = surface_laplacian(ucf,m)


V = TestFESpace(cube_model_2D, ReferenceFE(lagrangian,Float64,2); conformity=:H1)
U = TrialFESpace(V)

Jp(p) = x -> r1p_3D[p]⋅Jpanel1(x)
_Jcf = map(p->GenericField(Jp(p)),panel_ids)
Jcf = CellData.GenericCellField(_Jcf,Ω,PhysicalDomain())
Gcf = (Operation(transpose)(Jcf)) ⋅ Jcf
invG = Operation(inv)(Gcf)
measG = Operation((meas))(Gcf)
sqrtmeasG = Operation(sqrt)(measG)

stiffnes(u,v) = ∫( (Jcf ⋅(invG ⋅ ∇(v)))⋅ (Jcf ⋅invG ⋅ ∇(u) ) *(sqrtmeasG))dΩ
bound(v) = ∫( (v*(invG⋅ ∇(ucf))⋅n_Γ )*(sqrtmeasG) )dΓ

rhs = ucf + 1.0*(surface_laplacian(ucf,m))
sum(∫(rhs)dΩ)
helmholtz_biform(u,v) =  ∫( (u*v)* (sqrtmeasG) )dΩ - stiffnes(u,v)
helmholtz_liform(v) = ∫(  (rhs*v) * (sqrtmeasG) )dΩ - bound(v)

bb =  assemble_vector(bound,V)
∇(ucf)(get_cell_points(Γ))
n_Γ(get_cell_points(Γ))

op = AffineFEOperator(helmholtz_biform,helmholtz_liform,U,V)

uh = solve(LUSolver(),op)

# (uh-ucf)(pts)
e = l2(uh-ucf,dΩ)



u_p1 = (uh)(pts)[panel_ids.==1]
u_p2 = (ucf-uh)(pts)[panel_ids.==2]
u_p3 = (ucf-uh)(pts)[panel_ids.==3]
u_p4 = (ucf-uh)(pts)[panel_ids.==4]
u_p5 = (ucf-uh)(pts)[panel_ids.==5]
u_p6 = (ucf-uh)(pts)[panel_ids.==6]


A = get_matrix(op)
eigvals(Array(A))

quad = dΩ.quad
cell_map = get_cell_map(quad.trian)
cell_Jt = lazy_map(∇,cell_map)
cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)

meas(cell_Jtx[1][1])




function uαβ_scalar(p)
  function _h(αβ)
    if p == 2 || p == 6 || p == 4
      return -αβ[1]*αβ[2]
    else
      return αβ[1]*αβ[2]
    end
  end
end

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

ambient_model = get_ambient_model(manifold_model)
Ω_amb = Triangulation(ambient_model)

cell_field = map(p->GenericField(u_scalar_parametric2ambient(p,uαβ_scalar(p))),panel_ids)
ucf_ambient = CellData.GenericCellField(cell_field,Ω_amb,PhysicalDomain())
ucf_ambient(get_cell_points(Ω_amb))
writevtk(Ω_amb,dir*"/poisson_ambient",cellfields=["u"=>ucf_ambient],append=false)
