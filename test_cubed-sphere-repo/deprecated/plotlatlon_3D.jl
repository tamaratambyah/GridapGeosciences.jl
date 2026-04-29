using Gridap
using GridapGeosciences
using Gridap.Geometry
using DrWatson
using Plots

using MPI
using PartitionedArrays
using GridapGeosciences.Distributed
using GridapP4est
using GridapDistributed
using Test


include("convergence_tools.jl")
dir = datadir("plot_latlon")
!isdir(dir) && mkdir(dir)


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


num_horizontal_uniform_refinements = 5
num_vertical_uniform_refinements = 5
o3model = Parametric3DOctreeDistributedDiscreteModel(ranks,radius,thickness;
                                           num_horizontal_uniform_refinements=num_horizontal_uniform_refinements,
                                           num_vertical_uniform_refinements=num_vertical_uniform_refinements);


panel_model = o3model.parametric_dmodel


### Find the unique r values
_panel_model = panel_model.models.item
_panel_ids = get_panel_ids(_panel_model)
cmap = get_cell_map(get_grid(_panel_model))
ref_points = get_cell_ref_coordinates(_panel_model)
coords = lazy_map(evaluate,cmap,ref_points)
cell_geo_map = lazy_map(p -> ForwardMap(p), _panel_ids)
xyz = lazy_map(evaluate,cell_geo_map,coords)

fi = lazy_map(p->Cartesian2SphericalMap3D(),_panel_ids)
latlon_cell_geo_map = lazy_map(∘, fi, cell_geo_map)
θϕr = lazy_map(evaluate,latlon_cell_geo_map,coords)

Rs = []
for (i,θϕr) in enumerate(θϕr)
  r = map(x->x[3],θϕr)
  push!(Rs,r...)
end
unique_r = unique(round.(Rs,digits=8))

### test the unique r values equivalent to thicnkess/nc+1
n = _nc_vertical(panel_model)
_unique_r  = RADIUS .+ collect(range(0,THICKNESS,n+1))

@test all(unique_r .== _unique_r)


##### plot the 3D linear bousineq initial condition
a_e = 6.37e6/125 # m
d = 5000 #m
Lz = 20e3 #m
R = a_e # m radius
u_0 = 20 #m/s
Ωr = 7.292e-5 #1/s
c = 343 #m/s speed of sound
N = 0.01 #1/s bouyancy frequency
ztop = 10e3 #m
dΘ = 1 #K
TF = (3600*24)*10 # s

LH = a_e # m
LV = ztop/THICKNESS
c2 = c/LH # 1/s
τ = 1/c2#1/Ωr # s

_R = R/LH
_ztop = ztop/LV
_d = d/LV
_Lz = Lz/LV
_u_0 = u_0*τ/LH
_Ωr = Ωr*τ
_c = c*τ/LH
_N = N*τ
tF = TF/τ

function b0(xyz)
  x,y,z = xyz
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr

  θc = 2*π/3
  ϕc = 0.0

  k = sqrt(x^2 + y^2 + z^2) - _R

  r = _R*acos( sin(ϕc)*sin(ϕ) + cos(ϕc)*cos(ϕ)*cos(θ-θc)    )
  s = _d^2/(_d^2 + r^2)
  b = dΘ*s*sin( 2*π*k/_Lz  )
  b
end
b = panel_to_cartesian(b0)

Ω_panel = Triangulation(panel_model)
panel_ids = get_panel_ids(Ω_panel)
b_cf = ParametricCellField(b,Ω_panel,panel_ids)

owned_panel_ids = get_owned_panel_ids(panel_model)

latlon_cell_geo_map = map(owned_panel_ids) do pid
  cell_geo_map = lazy_map(p -> ForwardMap(p), pid)
  fi = lazy_map(p->Cartesian2SphericalMap3D(),pid)
  return lazy_map(∘, fi, cell_geo_map)
end

writevtk(Ω_panel,dir*"/3D_latlon_model",cellfields=["b"=>b_cf],append=false,geo_map=latlon_cell_geo_map,order=1)
