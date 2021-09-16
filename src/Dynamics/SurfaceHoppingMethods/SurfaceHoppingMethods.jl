
module SurfaceHoppingMethods

using DEDataArrays: DEDataVector
using StatsBase: sample, Weights

mutable struct SurfaceHoppingVariables{T,A,Axes,S} <: DEDataVector{T}
    x::ComponentVector{T,A,Axes}
    state::S
end

get_velocities(u::SurfaceHoppingVariables) = get_velocities(u.x)
get_positions(u::SurfaceHoppingVariables) = get_positions(u.x)
get_quantum_subsystem(u::SurfaceHoppingVariables) = get_quantum_subsystem(u.x)

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
abstract type SurfaceHopping <: Method end

function motion!(du, u, sim::AbstractSimulation{<:SurfaceHopping}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    dσ = get_quantum_subsystem(du)

    r = get_positions(u)
    v = get_velocities(u)
    σ = get_quantum_subsystem(u)

    set_state!(u, sim.method.state) # Make sure the state variables match, 

    velocity!(dr, v, r, sim, t)
    Calculators.update_electronics!(sim.calculator, r)
    acceleration!(dv, v, r, sim, t, sim.method.state)
    set_quantum_derivative!(dσ, v, σ, sim)
end

function check_hop!(u, t, integrator)::Bool
    sim = integrator.p
    evaluate_hopping_probability!(sim, u, get_proposed_dt(integrator))
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

const HoppingCallback = DiscreteCallback(check_hop!, execute_hop!; save_positions=(false, false))
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

function create_problem(u0, tspan, sim::AbstractSimulation{<:SurfaceHopping})
    set_state!(sim.method, u0.state)
    ODEProblem(motion!, u0, tspan, sim; callback=get_callbacks(sim))
end

include("fssh.jl")
include("iesh.jl")
include("rpsh.jl")

end # module