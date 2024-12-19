using Documenter
using DocumenterCitations
using NQCBase, NQCModels, NQCDistributions, NQCDynamics
using FrictionProviders

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

function find_all_files(directory)
    map(
        s -> joinpath(directory, s),
        sort(readdir(joinpath(@__DIR__, "src", directory)))
    )
end

@time makedocs(;
    plugins=[bib],
    sitename="NQCDynamics.jl",
    modules=[NQCDynamics, NQCDistributions, NQCModels, NQCBase, FrictionProviders],
    doctest=false,
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://nqcd.github.io/NQCDynamics.jl/stable/",
        assets=["assets/favicon.ico", "assets/citations.css"],
    ),
    authors="James Gardner and contributors.",
    pages=[
        "Introduction" => "index.md"
        "Getting started" => "getting_started.md"
        "Atoms" => "atoms.md"
        "Ensemble simulations" => "ensemble_simulations.md"
        "Saving and loading" => "saving_loading.md"
        "NQCModels.jl" => Any[
            "NQCModels/overview.md"
            "NQCModels/analyticmodels.md"
            "NQCModels/ase.md"
            "NQCModels/machinelearningmodels.md"
            "NQCModels/frictionmodels.md"
        ]
        "NQCDistributions.jl" => Any[
            "NQCDistributions/overview.md"
        ]
        "Initial conditions" => find_all_files("initialconditions")
        "Dynamics simulations" => Any[
            "dynamicssimulations/dynamicssimulations.md"
            find_all_files("dynamicssimulations/dynamicsmethods")
        ]
        "Outputs and Analysis" => find_all_files("output_and_analysis")
        "Examples" => find_all_files("examples")
        "integration_algorithms.md"
        "Developer documentation" => find_all_files("devdocs")
        "API" => Any[
            "NQCBase" => find_all_files("api/NQCBase")
            "NQCModels" => find_all_files("api/NQCModels")
            "NQCDistributions" => find_all_files("api/NQCDistributions")
            "NQCDynamics" => find_all_files("api/NQCDynamics")
        ]
        "References" => "references.md"
    ])


if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo="github.com/NQCD/NQCDynamics.jl",
        push_preview=true
    )
end
