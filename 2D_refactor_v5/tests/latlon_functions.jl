using Gridap
include("../src/initialise.jl")

struct SigmaMap2{A} <: Map
  r::A # sphere radius
end

""" map 3D point on sphere -> latlon (2D) """
function Gridap.Arrays.return_cache(k::SigmaMap2,cellx::AbstractArray{<:VectorValue{3,T}}) where {T}
  y = similar(cellx,VectorValue{2,T})
  return y #CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::SigmaMap2,cellx::AbstractArray{<:VectorValue{3}})
  # 3D point on sphere -> lat lon
  y = cache
  map!(x -> VectorValue( rem2pi(atan(x[2], x[1]),RoundNearest),
                                 atan(x[3], sqrt(x[1]*x[1] + x[2]*x[2]) )),
                        y, cellx)
  return y
end

function Gridap.Arrays.return_cache(k::SigmaMap2,x::VectorValue{3,T}) where {T}
  y = zero(VectorValue{2,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SigmaMap2,x::VectorValue{3})
  # 3D point on sphere -> lat lon
  # println("this func")
  y = cache
  y = VectorValue(  rem2pi(atan(x[2], x[1]),RoundNearest),
                           atan(x[3], sqrt(x[1]*x[1] + x[2]*x[2]) ))
  return y
end





manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

ambient_model = get_ambient_model(manifold_model)

ambient_grid = get_grid(ambient_model)
ambient_nodes = get_node_coordinates(ambient_model)
latlon_nodes = collect1d(lazy_map(SigmaMap2(r),ambient_nodes))

Ω_parametric = Triangulation(manifold_model)


latlon_grid = UnstructuredGrid(latlon_nodes,get_cell_node_ids(ambient_grid),
      get_reffes(ambient_grid),get_cell_type(ambient_grid),OrientationStyle(ambient_grid))

latlon_model = UnstructuredDiscreteModel(latlon_grid,get_grid_topology(manifold_model),get_face_labeling(manifold_model))

Ω_latlon = Triangulation(latlon_model)
pts_latlon = get_cell_points(Ω_latlon)

Ω_ambient = Triangulation(ambient_model)


function uθϕ(θϕ)
  θ,ϕ = θϕ
  cos(θ)
end

function u_latlon(p::Int)
  function _u(αβ)
    latlon_panel1 = evaluate(GnomonicMap(), αβ)
    sphere_panel1 = evaluate(Sigma(),latlon_panel1)
    sphere_panelp = evaluate( PanelRotationMap(rp1_3D[p]), sphere_panel1)
    latlon_panelp = evaluate(SigmaMap2(r),sphere_panelp)
    uθϕ(latlon_panelp)
  end
end

cell_field = map(p->GenericField(u_latlon(p)),panel_ids)
cf_parametric = CellData.GenericCellField(cell_field,Ω_parametric,PhysicalDomain())

dΩ = Measure(Ω_parametric,2)
H1 = FESpace(manifold_model,ReferenceFE(lagrangian,Float64,1), conformity=:H1)
biform(u,v) = ∫(u*v)dΩ
liform(v) = ∫(v*cf_parametric )dΩ

op = AffineFEOperator(biform,liform,H1,H1)
uh = solve(LUSolver(),op)

writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>uh],append=false)



pts_ambient = get_cell_points(Ω_ambient)

mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)
cf_mapped = lazy_map(Broadcasting(∘),get_data(uh),mapping)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)


cf_ambient(pts_ambient)
1;



















################################################################################
_perms = get_vertex_permutations(QUAD)
perms_1p = [_perms[1],_perms[1],_perms[1],_perms[4],[1,4,2,3],_perms[4]]
perms_p1 = [_perms[1],_perms[1],_perms[1],reverse(_perms[4]),reverse(_perms[5]),reverse(_perms[4])]


latlon_coords = get_cell_coordinates(latlon_model)
Panel1 = latlon_coords[panel_ids .==1]
Panel2 = latlon_coords[panel_ids .==2]
Panel3 = latlon_coords[panel_ids .==3]


# for cell in 1:24
cell = 24
  if Int(mod(cell,4)) > 0
    id = Int(mod(cell,4))
  else
    id = 4
  end
  panel1_latlon = Panel1[id]
  p = panel_ids[cell]
  ids = [p for i in 1:4]
  c = [0,0,-1,1]

  panel_latlon = lazy_map(AddtoLatLon(),panel1_latlon,ids,c)
  println(panel_latlon  ≈ latlon_coords[cell])
# end



struct AddtoLatLon <: Map
end


function Gridap.Arrays.return_cache(k::AddtoLatLon,cellθϕ1::VectorValue{2},panel_id::Int,c::Int)
  y = zero(VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::AddtoLatLon,cellθϕ1::VectorValue{2},panel_id::Int,c::Int)

  y = cache
  _R = [1 0 0
        0 0 -1
        0 1 0]
    R = TensorValue(_R) # rotation in axis by pi/2
    X1 = cos(cellθϕ1[1])*cos(cellθϕ1[2])
    Y1 = sin(cellθϕ1[1])*cos(cellθϕ1[2])
    Z1 = sin(cellθϕ1[2])
    XX1 = VectorValue(X1,Y1,Z1)

    XXR = R⋅XX1

    θr = rem2pi(atan(XXR[2], XXR[1]),RoundNearest)
    ϕr =  atan(XXR[3], sqrt(XXR[1]*XXR[1] + XXR[2]*XXR[2]) )


  if panel_id == 1
    y = VectorValue( cellθϕ1[1],cellθϕ1[2])
  elseif panel_id == 2
    # y = VectorValue( rem2pi( cellθϕ1[1] + c*π/2, RoundNearest ),abs(cellθϕ1[2]))
    y = VectorValue(cellθϕ1[1],cellθϕ1[2]+π/2)
  elseif panel_id == 3

    y = VectorValue( rem2pi( cellθϕ1[1] + π/2, RoundDown ),cellθϕ1[2])
  elseif panel_id == 4
    y = VectorValue( rem2pi(θr + π, RoundNearest ), ϕr)
  # elseif panel_id == 5
  #   y = VectorValue( rem2pi( cellθϕ1[1] + c*π, RoundDown ),-abs(cellθϕ1[2]))
  elseif panel_id == 6
    y = VectorValue( rem2pi( θr - π/2, RoundNearest ),ϕr)
  end
  return y
#
end



################################################################################
RT = FESpace(latlon_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
function uθϕ(θϕ)
  θ,ϕ = θϕ
  VectorValue(cos(θ),0.0)
end

uh = interpolate_everywhere(uθϕ,RT)
writevtk(Ω_latlon,dir*"/latlon_vector",cellfields=["u"=>uh],append=false)


################################################################################
## Girdalo map

shifts = [[0.0, 0.0 ], # front
          [0.0,π/2], # top
          [π/2,0.0], # left
          [π,0.0],  # back
          [0.0,-π/2], # bottom
          [3*π/2,0.0] # right
]

struct GirdaloMap{A} <: Map
  shifts::A
end


function Gridap.Arrays.return_cache(k::GirdaloMap,cellθϕ1::VectorValue{2,T},panel_id::Int) where {T}
  y = zero(VectorValue{2,T})
  z = 0.0
  w = 0.0
  return y,z,w
end

function Gridap.Arrays.evaluate!(cache,f::GirdaloMap,cellθϕ1::VectorValue{2},panel_id::Int)

  y,ϕ,θ = cache
  shift = f.shifts[panel_id]

  ϕ =  asin( sin(cellθϕ1[2])*cos(shift[2]) +  cos(cellθϕ1[1])*cos(cellθϕ1[2])*sin(shift[2]) )

  θ = rem2pi( shift[1] +
        atan( cos( ϕ)*sin(cellθϕ1[1]), ( cos(cellθϕ1[2])*cos(cellθϕ1[1])*cos(shift[2]) - sin(cellθϕ1[2])*sin(shift[2])  ) ),
                        RoundNearest)

  y = VectorValue(θ,ϕ)


  return y
#
end

latlon_coords = get_cell_coordinates(latlon_model)
Panel1 = latlon_coords[panel_ids .==1]

cell = 20
p = 5
p_ids = [p for i in 1:4]

panel1_latlon = Panel1[1]
panel_latlon = lazy_map(GirdaloMap(shifts),panel1_latlon,p_ids)

latlon_coords[cell]
