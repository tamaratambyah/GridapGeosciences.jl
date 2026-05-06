"""
This module tests the construction of GridapSolvers.MultilevelTools.ModelHierarchy
with parametric coarse models
"""


module HierarchyTest

using Gridap
using GridapGeosciences
using GridapSolvers
using Test

function main(model0,n_ref_lvls)
  mh = ModelHierarchy(model0,n_ref_lvls)
  @test isa(mh,ModelHierarchy)
  reffe = ReferenceFE(lagrangian,Float64,1)
  tests  = TestFESpace(mh,reffe)
  @test true
end


end ## module
