using Random
using LinearAlgebra
using ..Calculators: DiabaticFrictionCalculator

export MDEF
export TwoTemperatureMDEF

abstract type AbstractMDEF <: Method end

"""
$(TYPEDEF)

```math
dr = v dt\\\\
dv = -\\Delta U/M dt - \\Gamma v dt + \\sigma \\sqrt{2\\Gamma} dW
```
``\\Gamma`` is the friction tensor with units of inverse time.
For thermal dynamics we set ``\\sigma = \\sqrt{kT / M}``,
where ``T`` is the electronic temperature.

This is integrated using the BAOAB algorithm where the friction "O" step is performed
in the tensor's eigenbasis. See `src/dynamics/mdef_baoab.jl` for details.
"""
struct MDEF <: AbstractMDEF end

"""
$(TYPEDEF)

Same as standard MDEF but uses a function to determine the time-dependent temperature.
"""
struct TwoTemperatureMDEF <: AbstractMDEF
    temperature::Function
end

"""Gets the temperature as a function of time during MDEF."""
get_temperature(sim::Simulation{MDEF}, ::AbstractFloat) = sim.temperature
get_temperature(sim::Simulation{TwoTemperatureMDEF}, t::AbstractFloat) = sim.method.temperature(t)

"""
    acceleration!(dv, v, r, sim::Simulation{MDEF,<:DiabaticFrictionCalculator}, t)

Sets acceleration due to ground state force when using a `DiabaticFrictionModel`.
"""
function acceleration!(dv, v, r, sim::Simulation{MDEF,<:DiabaticFrictionCalculator}, t)
    Calculators.update_electronics!(sim.calculator, r)
    for i in axes(r, 2)
        for j in axes(r, 1)
            dv[j,i] = -sim.calculator.adiabatic_derivative[j,i][1,1] / sim.atoms.masses[i]
        end
    end
end

"""
    friction!(du, r, sim, t)

Evaluates friction tensor and provides variance of random force.
"""
function friction!(du, r, sim, t)
    Calculators.evaluate_friction!(sim.calculator, r)

    du.x[1] .= sim.calculator.friction
    for i in range(sim.atoms)
        for j in range(sim.atoms)
            for k=0:sim.DoFs-1
                du.x[1][i+k,j+k] /= sqrt(sim.atoms.masses[i] * sim.atoms.masses[j])
            end
        end
    end

    du.x[2] .= sqrt.(get_temperature(sim, t) ./ repeat(sim.atoms.masses; inner=sim.DoFs))
end

function create_problem(u0::ClassicalDynamicals, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    create_problem(u0.x, tspan, sim)
end

function create_problem(u0::ArrayPartition, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    DynamicalSDEProblem(acceleration!, velocity!, friction!, get_velocities(u0), get_positions(u0), tspan, sim)
end

select_algorithm(::AbstractSimulation{<:AbstractMDEF}) = MDEF_BAOAB()
