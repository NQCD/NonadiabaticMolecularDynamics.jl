using StructArrays: StructArray
using ComponentArrays: ComponentArrays
using LinearAlgebra: diagm, mul!
using NQCDynamics: RingPolymers
using NQCModels: nstates
using RingPolymerArrays: get_centroid
using NQCDistributions: ElectronicDistribution, Adiabatic, Diabatic, density_matrix
using ..Calculators: AbstractDiabaticCalculator, DiabaticCalculator, RingPolymerDiabaticCalculator

function set_quantum_derivative! end

function calculate_density_matrix_propagator!(sim::Simulation, v)
    V = sim.method.density_propagator

    V .= diagm(sim.calculator.eigen.values)
    for I in eachindex(v)
        @. V -= im * v[I] * sim.calculator.nonadiabatic_coupling[I]
    end
    return V
end

function calculate_density_matrix_propagator!(sim::RingPolymerSimulation, v)
    V = sim.method.density_propagator
    centroid_v = get_centroid(v)

    V .= diagm(sim.calculator.centroid_eigen.values)
    for I in eachindex(centroid_v)
        @. V -= im * centroid_v[I] * sim.calculator.centroid_nonadiabatic_coupling[I]
    end
    return V
end

function commutator!(out, A, B, tmp)
    mul!(out, A, B)
    mul!(tmp, B, A)
    out .-= tmp
    return nothing
end

get_quantum_subsystem(u::ComponentArrays.ComponentVector{T}) where {T} =
    StructArray{Complex{T}}((u.σreal, u.σimag))

function initialise_adiabatic_density_matrix(
    electronics::ElectronicDistribution{Diabatic},
    calculator::AbstractDiabaticCalculator,
    r
)

    diabatic_density = density_matrix(electronics, nstates(calculator))
    return transform_density!(diabatic_density, calculator, r, :to_adiabatic)
end

function initialise_adiabatic_density_matrix(
    electronics::ElectronicDistribution{Adiabatic},
    calculator::AbstractDiabaticCalculator,
    r
)

    adiabatic_density = density_matrix(electronics, nstates(calculator))
    return adiabatic_density
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
    eigs =  Calculators.get_eigen(calculator, r)
    return eigs.vectors
end

function evaluate_transformation(calculator::RingPolymerDiabaticCalculator, r)
    centroid_eigs = Calculators.get_centroid_eigen(calculator, r)
    return centroid_eigs.vectors
end
