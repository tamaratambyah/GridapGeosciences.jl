# module WaveEquationTestsSeq
include("../WaveEquation.jl")

## Serial model: 2D
n_ref_lvls = 4
models = get_refined_models(n_ref_lvls)
main(models)

# end
