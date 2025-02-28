using Documenter
using DocumenterCitations
using NQCBase, NQCModels, NQCDistributions, NQCDynamics
using FrictionProviders
using MACEModels

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

function find_all_files(directory)
    files = String[]
    for potential_file in sort(readdir(joinpath(@__DIR__, "src", directory)))
        if potential_file[1] != "." # Ensure we aren't adding any hidden files
            push!(files, joinpath(directory, potential_file))
        end
    end
    return files
end

@time makedocs(;
    plugins=[bib],
    sitename="NQCDynamics.jl",
    modules=[NQCDynamics, NQCDistributions, NQCModels, NQCBase, MACEModels, FrictionProviders],
    doctest=false,
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://nqcd.github.io/NQCDynamics.jl/stable/",
        assets=["assets/favicon.ico", "assets/citations.css"],
        size_threshold = 5000*1024,
        example_size_threshold = 8000*1024,
        size_threshold = 5000*1024,
        example_size_threshold = 8000*1024,
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
            "NQCModels/combining_models.md"
            "NQCModels/combining_models.md"
            "NQCModels/analyticmodels.md"
            "NQCModels/machinelearningmodels.md"
            "NQCModels/fullsizemodels.md"
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
