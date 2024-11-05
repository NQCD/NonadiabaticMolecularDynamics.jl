using StatsBase: mean
using NQCDynamics: masses
using NQCDistributions: ElectronicDistribution

export FSSH

"""
    FSSH{T} <: SurfaceHopping

Type for fewest-switches surface hopping

```jldoctest
Simulation{FSSH}(Atoms(:H), Free())

# output

Simulation{FSSH{Float64}}:
  Atoms{Float64}([:H], [1], [1837.4715941070515])
  Free(1)
```
"""
mutable struct FSSH{T} <: SurfaceHopping
    density_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    state::Int
    new_state::Int
    rescaling::Symbol
    tmp_complex_matrix::Matrix{Complex{T}}
    function FSSH{T}(states::Integer, rescaling::Symbol) where {T}
        density_propagator = zeros(states, states)
        hopping_probability = zeros(states)
        tmp_complex_matrix = zeros(Complex{T}, states, states)
        new{T}(density_propagator, hopping_probability, 0, 0, rescaling, tmp_complex_matrix)
    end
end

function Simulation{FSSH}(atoms::Atoms{T}, model::Model; rescaling=:standard, kwargs...) where {T}
    Simulation(atoms, model, FSSH{T}(NQCModels.nstates(model), rescaling); kwargs...)
end

function DynamicsMethods.DynamicsVariables(
    sim::AbstractSimulation{<:SurfaceHopping}, v, r, electronic::ElectronicDistribution
)
    σ = DynamicsUtils.initialise_adiabatic_density_matrix(electronic, sim.calculator, r)
    state = sample(Weights(diag(real.(σ))))
    return SurfaceHoppingVariables(ComponentVector(v=v, r=r, σreal=σ, σimag=zero(σ)), state)
end

function DynamicsUtils.acceleration!(dv, v, r, sim::AbstractSimulation{<:FSSH}, t, state)
    adiabatic_derivative = Calculators.get_adiabatic_derivative(sim.calculator, r)
    for I in eachindex(dv)
        dv[I] = -adiabatic_derivative[I][state, state]
    end
    DynamicsUtils.divide_by_mass!(dv, masses(sim))
    return nothing
end

"""
    evaluate_hopping_probability!(sim::Simulation{<:FSSH}, u, dt)

Evaluates the probability of hopping from the current state to all other states

# Implementation

- `σ` is Hermitan so the choice `σ[m,s]` or `σ[s,m]` is irrelevant; we take the real part.
- 'd' is skew-symmetric so here the indices are important.
"""
function evaluate_hopping_probability!(sim::AbstractSimulation{<:FSSH}, u, dt)
    v = DynamicsUtils.get_hopping_velocity(sim, DynamicsUtils.get_velocities(u))
    σ = DynamicsUtils.get_quantum_subsystem(u)
    s = sim.method.state
    r = DynamicsUtils.get_positions(u)
    d = DynamicsUtils.get_hopping_nonadiabatic_coupling(sim, r)

    fewest_switches_probability!(sim.method.hopping_probability, v, σ, s, d, dt)
end

function fewest_switches_probability!(probability, v, σ, s, d, dt)
    probability .= 0 # Set all entries to 0
    for m in axes(σ, 1)
        if m != s
            for I in eachindex(v)
                probability[m] += 2v[I]*real(σ[m,s]/σ[s,s])*d[I][s,m] * dt
            end
        end
    end

    clamp!(probability, 0, 1) # Restrict probabilities between 0 and 1
    cumsum!(probability, probability)
end

function select_new_state(sim::AbstractSimulation{<:FSSH}, u)

    random_number = rand()
    for (i, prob) in enumerate(sim.method.hopping_probability)
        if i != sim.method.state # Avoid self-hops
            if prob > random_number
                return i
            end
        end
    end
    return sim.method.state
end

"""
    unpack_states(sim)

Get the two states that we are hopping between.
"""
function unpack_states(sim::AbstractSimulation{<:FSSH})
    return (sim.method.new_state, sim.method.state)
end

function Estimators.diabatic_population(sim::AbstractSimulation{<:FSSH}, u)
    r = DynamicsUtils.get_positions(u)
    U = DynamicsUtils.evaluate_transformation(sim.calculator, r)

    σ = copy(DynamicsUtils.get_quantum_subsystem(u).re)
    σ[diagind(σ)] .= 0
    σ[u.state, u.state] = 1

    return diag(U * σ * U')
end

function Estimators.adiabatic_population(sim::AbstractSimulation{<:FSSH}, u)
    population = zeros(NQCModels.nstates(sim.calculator.model))
    population[u.state] = 1
    return population
end

function DynamicsUtils.classical_potential_energy(sim::Simulation{<:FSSH}, u)
    eigs = Calculators.get_eigen(sim.calculator, DynamicsUtils.get_positions(u))
    potential = eigs.values[u.state]
    return potential
end
