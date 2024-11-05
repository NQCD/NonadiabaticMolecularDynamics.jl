using LinearAlgebra: eigvecs, diag, diagind, dot
using StatsBase: sample, Weights
using RingPolymerArrays: get_centroid

using NQCDynamics.Calculators: RingPolymerDiabaticCalculator
using NQCDynamics: RingPolymerSimulation, RingPolymers

function RingPolymerSimulation{FSSH}(atoms::Atoms{T}, model::Model, n_beads::Integer; rescaling=:standard, kwargs...) where {T}
    RingPolymerSimulation(atoms, model, FSSH{T}(NQCModels.nstates(model), rescaling), n_beads; kwargs...)
end

function DynamicsMethods.motion!(du, u, sim::RingPolymerSimulation{<:SurfaceHopping}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    dσ = DynamicsUtils.get_quantum_subsystem(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    σ = DynamicsUtils.get_quantum_subsystem(u)

    set_state!(u, sim.method.state) # Make sure the state variables match, 

    DynamicsUtils.velocity!(dr, v, r, sim, t)
    Calculators.update_electronics!(sim.calculator, r)
    DynamicsUtils.acceleration!(dv, v, r, sim, t, sim.method.state)
    DynamicsUtils.apply_interbead_coupling!(dv, r, sim)
    DynamicsUtils.set_quantum_derivative!(dσ, u, sim)
end

function perform_rescaling!(
    sim::RingPolymerSimulation{<:SurfaceHopping}, velocity, γ, d
)
    for I in CartesianIndices(d)
        velocity[I,:] .-= γ * d[I] / masses(sim, I)
    end
    return nothing
end

function frustrated_hop_invert_velocity!(
    sim::RingPolymerSimulation{<:SurfaceHopping}, velocity, d
)
    #invert center of mass velocity of ring polymer
    dn = LinearAlgebra.normalize(d)
    vcom = sum(velocity, dims=3) / size(velocity,3)
    γ = dot(vcom[:,:,1],dn)
    for I in CartesianIndices(dn)
        velocity[I,:] .-= 2γ * dn[I]
    end
    return nothing
end

function DynamicsUtils.classical_potential_energy(sim::RingPolymerSimulation{<:FSSH}, u)
    all_eigs = Calculators.get_eigen(sim.calculator, DynamicsUtils.get_positions(u))
    potential = sum(eigs.values[u.state] for eigs in all_eigs)
    return potential
end

function DynamicsUtils.centroid_classical_potential_energy(sim::RingPolymerSimulation{<:FSSH}, u)
    centroid_eigs = Calculators.get_centroid_eigen(sim.calculator, DynamicsUtils.get_positions(u))
    potential = centroid_eigs.values[u.state]
    return potential
end
