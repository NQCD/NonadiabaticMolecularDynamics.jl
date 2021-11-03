using LinearAlgebra: diagind

using NonadiabaticModels: NonadiabaticModels
using NonadiabaticMolecularDynamics: Calculators
using .Calculators:
    AbstractDiabaticCalculator,
    DiabaticCalculator,
    RingPolymerDiabaticCalculator

abstract type StateType end
struct Diabatic <: StateType end
struct Adiabatic <: StateType end

abstract type ElectronicDistribution{S<:StateType} <: NonadiabaticDistribution end

struct SingleState{S} <: ElectronicDistribution{S}
    state::Int
    statetype::S
end
SingleState(state) = SingleState(state, Diabatic())

struct ElectronicPopulation{T,S} <: ElectronicDistribution{S}
    populations::Vector{T}
    statetype::S
end
ElectronicPopulation(state) = ElectronicPopulation(state, Diabatic())

function initialise_adiabatic_density_matrix(
    electronics::ElectronicDistribution{Diabatic},
    calculator::AbstractDiabaticCalculator,
    r
)

    diabatic_density = initialise_density_matrix(electronics, calculator)
    return transform_density!(diabatic_density, calculator, r, :to_adiabatic)
end

function initialise_adiabatic_density_matrix(
    electronics::ElectronicDistribution{Adiabatic},
    calculator::AbstractDiabaticCalculator,
    r
)

    adiabatic_density = initialise_density_matrix(electronics, calculator)
    return adiabatic_density
end

function initialise_density_matrix(
    electronics::ElectronicDistribution, calculator::AbstractDiabaticCalculator
)

    n = NonadiabaticModels.nstates(calculator.model)
    density = Matrix{eltype(calculator)}(undef, n, n)
    fill_density!(density, electronics)
    return density
end

function fill_density!(density::AbstractMatrix, electronics::SingleState)
    fill!(density, zero(eltype(density)))
    density[electronics.state, electronics.state] = 1
    return density
end

function fill_density!(density::AbstractMatrix, electronics::ElectronicPopulation)
    fill!(density, zero(eltype(density)))
    density[diagind(density)] .= electronics.populations
    return density
end

function transform_density!(
    density::AbstractMatrix, calculator::AbstractDiabaticCalculator, r, direction
)
    U = evaluate_transformation(calculator, r)
    if direction === :to_diabatic
        U = U'
    elseif !(direction === :to_adiabatic)
        throw(ArgumentError("`direction` $direction not recognised."))
    end
    density .= U' * density * U
    return density
end

function evaluate_transformation(calculator::DiabaticCalculator, r)
    Calculators.evaluate_potential!(calculator, r)
    Calculators.eigen!(calculator)
    return calculator.eigenvectors
end

function evaluate_transformation(calculator::RingPolymerDiabaticCalculator, r)
    Calculators.evaluate_centroid_potential!(calculator, r)
    Calculators.centroid_eigen!(calculator)
    return calculator.centroid_eigenvectors
end
