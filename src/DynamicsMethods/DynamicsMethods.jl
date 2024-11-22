"""
This module contains functions and types necessary for performing
nonadiabatic molecular dynamics.

Dynamics is performed using [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/).
As such, this module is centered around the implementation of the functions
necessary to integrate the dynamics.

For deterministic Hamiltonian methods, the central function is [`DynamicsMethods.motion!`](@ref),
which is the inplace form of the function
to be integrated by [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/).

Further, methods that have discontinuities, such as surface hopping, use the
[callback interface](https://diffeq.sciml.ai/stable/features/callback_functions/#callbacks)
provided by `DifferentialEquations.jl`.
"""
module DynamicsMethods

using DiffEqBase: DiffEqBase
using UnitfulAtomic: austrip
using Reexport: @reexport
using OrdinaryDiffEq: OrdinaryDiffEq
using ComponentArrays: ComponentVector

using NQCDynamics: AbstractSimulation, Simulation, RingPolymerSimulation
using NQCBase: austrip_kwargs

export DynamicsVariables
export run_trajectory

"""
Each type of dynamics subtypes `Method` which is passed to
the `AbstractSimulation` as a parameter to determine the type of
dynamics desired.

julia> ttree(DynamicsMethods.Method)
"""
abstract type Method end

"""
    motion!(du, u, sim, t)
    
As per `DifferentialEquations.jl`, this function is implemented for
each method and defines the time derivatives of the `DynamicalVariables`.

We require that each implementation ensures `du` and `u` are subtypes
of `DynamicalVariables` and `sim` subtypes `AbstractSimulation`.
"""
function motion! end

"Provides the DEProblem for each type of simulation."
create_problem(u0, tspan, sim) =
    OrdinaryDiffEq.ODEProblem(motion!, u0, tspan, sim; callback=get_callbacks(sim))

"Choose a default algorithm for solving the differential equation."
select_algorithm(::AbstractSimulation) = OrdinaryDiffEq.VCABM5()

"Select the default callbacks for this simulation type."
get_callbacks(::AbstractSimulation) = nothing

"""
    DynamicsVariables(::AbstractSimulation, args...)

For each dynamics method this function is implemented to provide the variables for the
dynamics in the appropriate format.

By default, `DynamicsVariables` is set up for the classical case and takes
`sim`, `v`, `r` as arguments and returns a `ComponentVector(v=v, r=r)`
which is used as a container for the velocities and positions during
classical dynamics.
"""
DynamicsVariables(::AbstractSimulation, v, r) = ComponentVector(v=v, r=r)

run_trajectory(args...; kwargs...) = throw(error("""
`run_trajectory` has been replaced, use `run_dynamics` instead.
`run_dynamics` unifies single trajectory and ensemble simulations into one function.
Refer to the `run_dynamics` docstring for more information.
"""
))

include("electronic_dynamics.jl")

include("ClassicalMethods/ClassicalMethods.jl")
@reexport using .ClassicalMethods

include("SurfaceHoppingMethods/SurfaceHoppingMethods.jl")
@reexport using .SurfaceHoppingMethods

include("MappingVariableMethods/MappingVariableMethods.jl")
@reexport using .MappingVariableMethods

include("EhrenfestMethods/EhrenfestMethods.jl")
@reexport using .EhrenfestMethods

include("IntegrationAlgorithms/IntegrationAlgorithms.jl")

end # module