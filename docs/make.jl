using Documenter

# Examples
include("examples.jl")

# Changelog
# cp(string(@__DIR__,"/../NEWS.md"), string(@__DIR__,"/src/changelog.md"))


# Make Docs
DocMeta.setdocmeta!(GridapGeosciences, :DocTestSetup, :(using GridapGeosciences); recursive=true)

makedocs(;
    modules = [GridapGeosciences],
    authors = "Tamara A. Tambyah <tamara.tambyah@monash.edu>, Alberto F. Martin <alberto.f.martin@anu.edu.au>, Santiago Badia <santiago.badia@monash.edu>, David Lee <david.lee@bom.gov.au>",
    # repo = "https://github.com/gridap/GridapGeosciences.jl/blob/{commit}{path}#{line}",
    repo = "https://github.com/tamaratambyah/GridapGeosciences.jl/tree/2D_refactor_v6",
    sitename = "GridapGeosciences.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://gridap.github.io/GridapGeosciences.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        # "Modules" => [
        #     "SolverInterfaces" => "SolverInterfaces.md",
        #     "MultilevelTools" => "MultilevelTools.md",
        #     "LinearSolvers" => "LinearSolvers.md",
        #     "NonlinearSolvers" => "NonlinearSolvers.md",
        #     "BlockSolvers" => "BlockSolvers.md",
        #     "PatchBasedSmoothers" => "PatchBasedSmoothers.md",
        # ],
        # "Extensions" => [
        #     "GridapP4est.jl" => "Extensions/GridapP4est.md",
        #     "GridapPETSc.jl" => "Extensions/GridapPETSc.md",
        #     "IterativeSolvers.jl" => "Extensions/IterativeSolvers.md",
        #     "Pardiso.jl" => "Extensions/Pardiso.md",
        # ],
        "Examples" => [
            "Laplace-Beltrami" => "Examples/LaplaceBeltrami.md",
            # "Navier-Stokes" => "Examples/NavierStokes.md",
            # "Stokes (GMG)" => "Examples/StokesGMG.md",
            # "Navier-Stokes (GMG)" => "Examples/NavierStokesGMG.md",
            # "Darcy (GMG)" => "Examples/DarcyGMG.md",
        ],
        # "Changelog" => "changelog.md",
    ],
    warnonly = [:doctest,:example_block,:eval_block],
    clean = true,
)

deploydocs(;
    # repo="github.com/gridap/GridapGeosciences.jl",
    repo="https://github.com/tamaratambyah/GridapGeosciences.jl",
    devbranch="2D_refactor_v6",
)
