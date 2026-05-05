# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

This is the first public release of GridapGeosciences. Non-exhaustive list of new features:

**Cubed sphere discrete models**
- `CubedSphere2DParametricDiscreteModel` is a serial implementation of the two dimensional cubed sphere, parametric model
- `CubedSphere2DParametricDistributedDiscreteModel` is a distributed implementation of the two dimensional cubed sphere, parametric model
- `CubedSphere2DParametricOctreeDistributedDiscreteModel` is a distributed, adaptive implementation of the two dimensional cubed sphere, parametric model, provided by `p4est` through `GridapP4est.jl`
- `CubedSphere3DParametricOctreeDistributedDiscreteModel` is a distributed, adaptive implementation of the three dimensional cubed sphere, parametric model, provided by `p4est` through `GridapP4est.jl`

**Continuous vector-valued Lagrangian finite elements**
- We provide a serial and distributed implementation of vector-valued Lagrangrian finite elements on manifold.

**Time integrator for differential algrabic equations**
- `DAEFEOperator` extends the ODE framework in `Gridap` to differential algebraic equations that arise in atmospheric systems like the shallow water equations. Currently implemented for explicit Runge Kutta methods

**Mapped vtk files**
- `writevtk_with_cell_geomap(geo_map::AbstractArray,...)` and `createvtk_with_cell_geomap(geo_map::AbstractArray,...)` extend the `writevtk(...)` functionality of `Gridap` by applying a cell-wise geometrical map the triangulation for visualisation purposes only. The evaluation of the fields remains on the input triangulation.
 


## [0.4.0] - 2023-10-17

A changelog is not maintained for this version.
Bump dependencies to the latest version


## [0.1.0] - 2021-12-09

A changelog is not maintained for this version.
Tagging the repo right before implementing support for distributed computing

