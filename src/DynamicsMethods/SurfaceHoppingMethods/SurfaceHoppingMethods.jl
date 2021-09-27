
module SurfaceHoppingMethods

using DEDataArrays: DEDataArrays
using ComponentArrays: ComponentVector
using DiffEqBase: DiffEqBase
using LinearAlgebra: lmul!
using OrdinaryDiffEq: OrdinaryDiffEq

using NonadiabaticMolecularDynamics:
    NonadiabaticMolecularDynamics,
    AbstractSimulation,
    Simulation,
    Calculators,
    DynamicsMethods,
    DynamicsUtils,
    Estimators,
    ndofs
using NonadiabaticModels: NonadiabaticModels, Model
using NonadiabaticDynamicsBase: Atoms

mutable struct SurfaceHoppingVariables{T,A,Axes,S} <: DEDataArrays.DEDataVector{T}
    x::ComponentVector{T,A,Axes}
    state::S
end

DynamicsUtils.get_velocities(u::SurfaceHoppingVariables) = DynamicsUtils.get_velocities(u.x)
DynamicsUtils.get_positions(u::SurfaceHoppingVariables) = DynamicsUtils.get_positions(u.x)
DynamicsUtils.get_quantum_subsystem(u::SurfaceHoppingVariables) = DynamicsUtils.get_quantum_subsystem(u.x)

"""
Abstract type for all surface hopping methods.

Surface hopping methods follow the structure set out in this file.
The nuclear and electronic variables are propagated by the `motion!` function.
The surface hopping procedure is handled by the `HoppingCallback` which
uses the functions `check_hop!` and `execute_hop!` as its `condition` and `affect!`.

To add a new surface hopping scheme, you must create a new struct
and define methods for `evaluate_hopping_probability!`, `select_new_state`,
and `rescale_velocity!`.

See `fssh.jl` for an example implementation.
"""
abstract type SurfaceHopping <: DynamicsMethods.Method end

function DynamicsMethods.motion!(du, u, sim::AbstractSimulation{<:SurfaceHopping}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    dσ = DynamicsUtils.get_quantum_subsystem(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    σ = DynamicsUtils.get_quantum_subsystem(u)

    set_state!(u, sim.method.state) # Make sure the state variables match, 

    DynamicsUtils.velocity!(dr, v, r, sim, t)
    Calculators.update_electronics!(sim.calculator, r)
    acceleration!(dv, v, r, sim, t, sim.method.state)
    set_quantum_derivative!(dσ, v, σ, sim)
end

function set_quantum_derivative!(dσ, v, σ, sim::AbstractSimulation{<:SurfaceHopping})
    V = DynamicsUtils.calculate_density_matrix_propagator!(sim, v)
    DynamicsUtils.commutator!(dσ, V, σ, sim.calculator.tmp_mat_complex1)
    lmul!(-im, dσ)
end

function check_hop!(u, t, integrator)::Bool
    sim = integrator.p
    evaluate_hopping_probability!(sim, u, OrdinaryDiffEq.get_proposed_dt(integrator))
    set_new_state!(sim.method, select_new_state(sim, u))
    return sim.method.new_state != sim.method.state
end

function execute_hop!(integrator)
    sim = integrator.p
    if rescale_velocity!(sim, integrator.u)
        set_state!(integrator.u, sim.method.new_state)
        set_state!(sim.method, sim.method.new_state)
    end
    return nothing
end

set_state!(container, new_state::Integer) = container.state = new_state
set_state!(container, new_state::AbstractVector) = copyto!(container.state, new_state)
set_new_state!(container, new_state::Integer) = container.new_state = new_state
set_new_state!(container, new_state::AbstractVector) = copyto!(container.new_state, new_state)

const HoppingCallback = DiffEqBase.DiscreteCallback(check_hop!, execute_hop!;
                                                    save_positions=(false, false))

get_callbacks(::AbstractSimulation{<:SurfaceHopping}) = HoppingCallback

"""
This function should set the field `sim.method.hopping_probability`.
"""
function evaluate_hopping_probability!(::AbstractSimulation{<:SurfaceHopping}, u, dt) end

"""
This function should return the desired state determined by the probability.
Should return the original state if no hop is desired.
"""
function select_new_state(::AbstractSimulation{<:SurfaceHopping}, u) end

"""
This function should modify the velocity and return a `Bool` that determines
whether the state change should take place.

This only needs to be implemented if the velocity should be modified during a hop.
"""
rescale_velocity!(::AbstractSimulation{<:SurfaceHopping}, u) = true

function DynamicsMethods.create_problem(u0, tspan, sim::AbstractSimulation{<:SurfaceHopping})
    set_state!(sim.method, u0.state)
    OrdinaryDiffEq.ODEProblem(DynamicsMethods.motion!, u0, tspan, sim; callback=get_callbacks(sim))
end

include("fssh.jl")
include("iesh.jl")
include("rpsh.jl")

end # module