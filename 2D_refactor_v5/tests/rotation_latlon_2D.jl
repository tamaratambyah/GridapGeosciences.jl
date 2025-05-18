using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using LinearAlgebra, FillArrays, StaticArrays

include("../src/initialise.jl")

_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
model = Adaptivity.refine(_model)
panel_ids = get_panel_ids(model)

X_cell_coords = get_cell_coordinates(get_ambient_model(model))

parametric_cell_coords =  get_grid(model).parametric_cell_coords

θϕ_cell_coords = lazy_map(SigmaMap(r),X_cell_coords)



shifts = [[0.0, 0.0 ],
          [0.0,π/2],
          [π/2,0.0],
          [π,0.0], #[π,0.0],
          [0.0,-π/2],
          [3*π/2,0.0]
]

_perms = get_vertex_permutations(QUAD)
perms_1p = [_perms[1],_perms[1],_perms[1],_perms[4],_perms[5],_perms[4]]
perms_p1 = [_perms[1],_perms[1],_perms[1],reverse(_perms[4]),reverse(_perms[5]),reverse(_perms[4])]


p = 4

shift = shifts[p]

panel1 = θϕ_cell_coords[panel_ids.==1]

ϕs = map(θϕ_cell_coords[1]) do x
  asin( sin(x[2])*cos(shift[2]) +  cos(x[1])*cos(x[2])*sin(shift[2]) )
end

θs = map(θϕ_cell_coords[1],ϕs) do x,phi
  rem2pi( shift[1] +
          atan( cos(phi)*sin(x[1]), ( cos(x[2])*cos(x[1])*cos(shift[2]) - sin(x[2])*sin(shift[2])  ) ),
  RoundDown)
end

θϕ = map(θs,ϕs) do θ,ϕ
  VectorValue(θ,ϕ)
end
θϕ_cell_coords[panel_ids.== p][1]


θϕ[perms_1p[p]] ≈ θϕ_cell_coords[panel_ids.== p][1]
θϕ_cell_coords[p][perms_p1[p]] ≈ θϕ

struct ParametricToLatLonMap{A,B} <: Map
  shifts::A
  perms::B
end


function Gridap.Arrays.return_cache(k::ParametricToLatLonMap,cellθϕ1::AbstractArray{<:VectorValue{2,T}},panel_id::Int) where {T}
  y = similar(cellθϕ1)
  z = similar(cellθϕ1,T)
  w = similar(cellθϕ1,T)
  return y,z,w
end

function Gridap.Arrays.evaluate!(cache,f::ParametricToLatLonMap,cellθϕ1::AbstractArray{<:VectorValue{2}},panel_id::Int)

  y,ϕs,θs = cache
  shift = f.shifts[panel_id]
  perm = f.perms[panel_id]

  map!(x-> asin( sin(x[2])*cos(shift[2]) +  cos(x[1])*cos(x[2])*sin(shift[2]) ),
    ϕs,cellθϕ1)

  map!((x,phi)->rem2pi( shift[1] +
                        atan( cos(phi)*sin(x[1]), ( cos(x[2])*cos(x[1])*cos(shift[2]) - sin(x[2])*sin(shift[2])  ) ),
                        RoundDown),
      θs,cellθϕ1,ϕs)

  map!( (i,j)-> VectorValue(i,j), y, θs,ϕs)

  # return y[perm]
  return y
#
end

θϕ1 =  lazy_map(GnomonicMap(),parametric_cell_coords)

z = lazy_map(ParametricToLatLonMap(shifts,perms_1p),θϕ1,panel_ids)
_X_cell_coords = lazy_map(SigmaMap(r),z)

# cache = array_cache(z)
# bm() = lazy_collect(cache,z)
# @benchmark bm()

(z .≈ θϕ_cell_coords)

cell = 24
z[cell]
θϕ_cell_coords[cell]



X_cell_coords[cell]
_X_cell_coords[cell]
