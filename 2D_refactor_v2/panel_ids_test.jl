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
using BenchmarkTools

include("initialise.jl")


ids = get_panel_ids(model)
ref_ids = get_panel_ids(ref_model)
ref_ref_ids = get_panel_ids(ref_ref_model)
ref_ref_ref_ids = get_panel_ids(ref_ref_ref_model)

_ref_ref_ref_ids = _get_panel_ids(ref_ref_ref_model)

bm1() = get_panel_ids(ref_ref_ref_model)
@benchmark bm1()

bm2() = _get_panel_ids(ref_ref_ref_model)
@benchmark bm2()


@allocated get_panel_ids(ref_ref_ref_model)
@allocated _get_panel_ids(ref_ref_ref_model)
