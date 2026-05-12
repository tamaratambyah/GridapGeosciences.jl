include("../AmbientTransientWaveEquation.jl")

## Serial model: 2D
n_ref_lvls = 3
radius = 1.0
models = get_ambient_refined_models(n_ref_lvls,radius)
AmbientTransientWaveEquationTests.main(models[1])
