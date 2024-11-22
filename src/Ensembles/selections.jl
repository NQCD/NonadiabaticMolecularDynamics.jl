
using NQCDistributions: DynamicalDistribution, ProductDistribution

"""
Check the tree structure of the type.
julia> ttree(Ensembles.AbstractSelection)
"""

abstract type AbstractSelection end

"""
Select the initial conditions from the distribution in order. 
"""
struct OrderedSelection{D,I} <: AbstractSelection
    "Distribution that is sampled."
    distribution::D
    indices::I
end

"""
Obtain initial conditions by randomly sampling the distribution.
"""
struct RandomSelection{D} <: AbstractSelection
    "Distribution that is sampled."
    distribution::D
end

function Selection(distribution, selection::AbstractVector, trajectories)
    if trajectories > 1
        @info "Sampling the provided distribution in the range $selection." 
    end
    OrderedSelection(distribution, selection)
end

function Selection(distribution, ::Any, trajectories)
    if trajectories > 1
        @info "Sampling randomly from provided distribution." 
    end
    RandomSelection(distribution)
end

function (select::OrderedSelection)(prob, i, repeat)
    j = select.indices[i]
    u0 = sample_distribution(prob.p, select.distribution, j)
    DynamicsMethods.create_problem(u0, prob.tspan, prob.p)
end

function sample_distribution(sim::AbstractSimulation, distribution::DynamicalDistribution, i)
    u = distribution[i]
    DynamicsMethods.DynamicsVariables(sim, u.v, u.r)
end

function sample_distribution(sim::AbstractSimulation, distribution::ProductDistribution, i)
    u = distribution.nuclear[i]
    DynamicsMethods.DynamicsVariables(sim, u.v, u.r, distribution.electronic)
end

function (select::RandomSelection)(prob, i, repeat)
    u0 = sample_distribution(prob.p, select.distribution)
    DynamicsMethods.create_problem(u0, prob.tspan, prob.p)
end

function sample_distribution(sim::AbstractSimulation, distribution::DynamicalDistribution)
    u = rand(distribution)
    DynamicsMethods.DynamicsVariables(sim, u.v, u.r)
end

function sample_distribution(sim::AbstractSimulation, distribution::ProductDistribution)
    u = rand(distribution.nuclear)
    DynamicsMethods.DynamicsVariables(sim, u.v, u.r, distribution.electronic)
end

sample_distribution(::AbstractSimulation, distribution) = copy(distribution)

