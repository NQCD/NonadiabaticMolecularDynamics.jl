using Documenter
using DocumenterCitations
using NQCBase, NQCModels, NQCDynamics
using CubeLDFAModel, NNInterfaces

DocMeta.setdocmeta!(NQCDynamics, :DocTestSetup, :(using NQCDynamics); recursive=true)
DocMeta.setdocmeta!(NQCModels, :DocTestSetup, :(using NQCModels, Symbolics); recursive=true)
DocMeta.setdocmeta!(NQCBase, :DocTestSetup, :(using NQCBase); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), sorting=:nyt)

function find_all_files(directory)
    map(
        s -> joinpath(directory, s),
        sort(readdir(joinpath(@__DIR__, "src", directory)))
    )
end

@time makedocs(
    bib,
    sitename = "NQCDynamics.jl",
    modules = [NQCDynamics, NQCModels, NQCBase, CubeLDFAModel],
    strict = false,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://nqcd.github.io/NQCDynamics.jl/stable/",
        assets = ["assets/favicon.ico"],
        ansicolor = true,
        ),
    authors = "James Gardner and contributors.",
    pages = [
        "Introduction" => "index.md"
        "Getting started" => "getting_started.md"
        "Atoms" => "atoms.md"
        "Ensemble simulations" => "ensemble_simulations.md"
        "Saving and loading" => "saving_loading.md"
        "NQCModels.jl" => Any[
            "NQCModels/overview.md"
            "NQCModels/analyticmodels.md"
            "NQCModels/ase.md"
            "NQCModels/neuralnetworkmodels.md"
            "NQCModels/frictionmodels.md"
        ]
        "Initial conditions" => Any[
            "initialconditions/dynamicaldistribution.md"
            find_all_files("initialconditions/samplingmethods")
        ]
        "Dynamics simulations" => Any[
            "dynamicssimulations/dynamicssimulations.md"
            find_all_files("dynamicssimulations/dynamicsmethods")
        ]
        "Examples" => find_all_files("examples")
        "Developer documentation" => find_all_files("devdocs")
        "API" => Any[
            "NQCBase" => find_all_files("api/NQCBase")
            "NQCModels" => find_all_files("api/NQCModels")
            "NQCDynamics" => find_all_files("api/NQCDynamics")
        ]
        "References" => "references.md"
    ])


if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/NQCD/NQCDynamics.jl",
        push_preview=true
    )
end
