using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Helpers
using Gridap.CellData
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools
include("helpers.jl")
include("maps/metric_maps.jl")

_gx = fill( [TensorValue{2,2}(1.0,2.0,3.0,4.0) for i in 1:10], 10)
_bx = fill( [1.0 for i in 1:10], 10)

# meass = lazy_map(MetricMeasure(),_gx)
# cache = array_cache(meass)
# bm() = lazy_collect(cache,meass)
# @benchmark bm()

# invv = lazy_map(MetricInverse(),_gx)
# cache = array_cache(invv)
# bm1() = lazy_collect(cache,invv)
# @benchmark bm1()

# out = lazy_map(LazyMult(), _bx, meass)
# cache = array_cache(out)
# bm1() = lazy_collect(cache,out)
# @benchmark bm1()

########################

metric_func(x) = TensorValue{2,2}(1.0, x[2]^2, (tan(x[1]))^2, (cos(x[1]))^2 )
sq_meas_func(x) = sqrt(meas(metric_func(x)))
inv_meas_func(x) = inv(metric_func(x))

parametric_model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π,π,-π/2,π/2),(8,8)))
Ω = Triangulation(parametric_model)

ref_pts = get_cell_points(Ω)
phys_pts = CellPoint(get_array(ref_pts),Ω,PhysicalDomain())
DomainStyle(phys_pts)
DomainStyle(ref_pts)

metric = CellField(metric_func,Ω)
DomainStyle(metric)
sq_meas = CellField(sq_meas_func,Ω)
inv_metric = CellField(inv_meas_func,Ω)


a = CellField(x->x[1],Ω)
b = CellField(x->VectorValue(x[1],x[2]),Ω)

surf_grad = inv_metric⋅ gradient(a)
surf_grad(phys_pts)./1

f = sq_meas*b
f(phys_pts)

surf_div = 1/sq_meas * divergence(f)
surf_div(phys_pts)

########################### debuggin


# op_sq_meas = Operation(x->sqrt(meas(x)))
# sq_meas = Gridap.CellData.OperationCellField(op_sq_meas,metric)
# sq_meas = op_sq_meas(metric)




# inv_metric = Gridap.CellData.OperationCellField(Operation(inv),metric)


# g_phys_cf = metric
# g_ref_cf = change_domain(metric,PhysicalDomain(),ReferenceDomain())
# DomainStyle(g_phys_cf)
# DomainStyle(g_ref_cf)

# gx_phys = g_phys_cf(phys_pts)
# gx_ref = g_ref_cf(ref_pts)
# test_coords(gx_phys,gx_ref)

# meass_phys = sq_meas(phys_pts)
# meass_ref = sq_meas(ref_pts)
# test_coords(meass_phys,meass_ref)

# inv_phys = lazy_map(MetricInverse(),gx_phys)
# inv_ref = lazy_map(MetricInverse(),gx_ref)
# test_coords(inv_phys,inv_ref)

# ginv_phys_cf = CellData.similar_cell_field(g_phys_cf,inv_phys,Ω,PhysicalDomain())
# ginv_ref_cf = CellData.similar_cell_field(g_ref_cf,inv_ref,Ω,ReferenceDomain())
