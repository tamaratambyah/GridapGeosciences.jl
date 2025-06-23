using Gridap
include("../src/initialise.jl")


# manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(1),cube)
manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
ref_manifold_model = Adaptivity.refine(manifold_model)
ref_ref_manifold_model = Adaptivity.refine(ref_manifold_model)

num_point_dims(manifold_model) == num_point_dims(ref_manifold_model) == num_point_dims(ref_ref_manifold_model)

@test is_child(ref_ref_manifold_model,ref_manifold_model)
@test is_child(ref_manifold_model,manifold_model)

# models = [ref_ref_manifold_model, ref_manifold_model, manifold_model]
# mh = ModelHierarchy(models)

writevtk(get_ambient_grid(get_grid(manifold_model)),dir*"/grid",append=false)
writevtk(get_ambient_grid(get_grid(ref_manifold_model)),dir*"/ref_grid",append=false)
writevtk(get_ambient_grid(get_grid(ref_ref_manifold_model)),dir*"/ref_ref_grid",append=false)

Ω = Triangulation(manifold_model)
ref_Ω = Triangulation(ref_manifold_model)
ref_ref_Ω = Triangulation(ref_ref_manifold_model)
