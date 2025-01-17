using DrWatson
using Gridap
using Gridap.Geometry
mapp(x) = Point(2*x[1],4*x[2],-1*x[3])
model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
# grid = get_grid(model)

_model = UnstructuredDiscreteModel(cube_grid)

Vh = FESpace(_model,
            ReferenceFE(lagrangian,VectorValue{3,Float64},1),
            conformity=:H1)
FE_map = interpolate(mapp,Vh)

ref_model = Gridap.Adaptivity.refine(_model)
ref_Vh = FESpace(ref_model,
            ReferenceFE(lagrangian,VectorValue{3,Float64},1),
            conformity=:H1)
ref_FE_map = interpolate(FE_map,ref_Vh)
