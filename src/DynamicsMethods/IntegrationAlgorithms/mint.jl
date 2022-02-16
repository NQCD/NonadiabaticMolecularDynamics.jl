
using UnPack: @unpack
using MuladdMacro: @muladd
using StaticArrays: SMatrix
using OrdinaryDiffEq: OrdinaryDiffEq
using LinearAlgebra: Hermitian, tr, Eigen, dot

using NQCDynamics.DynamicsMethods: MappingVariableMethods
using NQCModels: nstates, NQCModels

"""
    MInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm

Second order symplectic momentum integral algorithm.

# Reference

[J. Chem. Phys. 148, 102326 (2018)](https://doi.org/10.1063/1.5005557)
"""
struct MInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

OrdinaryDiffEq.isfsal(::MInt) = false

mutable struct MIntCache{uType,T} <: OrdinaryDiffEq.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    C::Matrix{T}
    D::Matrix{T}
    Γ::Matrix{T}
    Ξ::Matrix{T}
    tmp_mat::Matrix{T}
    tmp_vec1::Vector{T}
    tmp_vec2::Vector{T}
end

function OrdinaryDiffEq.alg_cache(::MInt,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    n = NQCModels.nstates(p.calculator.model)
    C = zeros(n,n)
    D = zeros(n,n)
    Γ = zeros(n,n)
    Ξ = zeros(n,n)
    tmp_mat = zeros(n,n)
    tmp_vec1 = zeros(n)
    tmp_vec2 = zeros(n)
    MIntCache(u, uprev, tmp, C, D, Γ, Ξ, tmp_mat, tmp_vec1, tmp_vec2)
end

function OrdinaryDiffEq.initialize!(_, ::MIntCache) end

@muladd function OrdinaryDiffEq.perform_step!(integrator, cache::MIntCache, repeat_step=false)
    @unpack dt,uprev,u,p = integrator
    @unpack tmp, Γ, Ξ = cache

    calc = p.calculator

    rtmp = DynamicsUtils.get_positions(tmp)
    vprev = DynamicsUtils.get_velocities(uprev)
    rprev = DynamicsUtils.get_positions(uprev)
    v = DynamicsUtils.get_velocities(u)
    r = DynamicsUtils.get_positions(u)
    X = MappingVariableMethods.get_mapping_positions(u)
    P = MappingVariableMethods.get_mapping_momenta(u)

    copyto!(v, vprev)
    copyto!(r, rprev)

    step_A!(rtmp, rprev, dt/2, v)

    propagate_mapping_variables!(cache, calc, X, P, rtmp, dt)

    eigs = Calculators.get_eigen(calc, rtmp)
    adiabatic_derivative = Calculators.get_adiabatic_derivative(calc, rtmp)
    state_independent_derivative = zero(rtmp)
    NQCModels.state_independent_derivative!(calc.model, state_independent_derivative, rtmp)
    X = MappingVariableMethods.get_mapping_positions(u)
    P = MappingVariableMethods.get_mapping_momenta(u)
    ∂V = Calculators.get_derivative(calc, rtmp)
    for i=1:natoms(p)
        for j=1:ndofs(p)

            set_gamma!(Γ, adiabatic_derivative[j,i], eigs.values, dt)
            transform_matrix!(Γ, eigs.vectors, cache.tmp_mat)
            E = Γ

            set_xi!(Ξ, adiabatic_derivative[j,i], eigs.values, dt)
            transform_matrix!(Ξ, eigs.vectors, cache.tmp_mat)
            F = Ξ

            force = get_mapping_nuclear_force(X, P, E, F, cache.tmp_vec1)
            v[j,i] -= force / p.atoms.masses[i]
            v[j,i] -= state_independent_derivative[j,i] / p.atoms.masses[i] * dt
            v[j,i] += tr(∂V[j,i]) * p.method.γ / p.atoms.masses[i] * dt / 2
        end
    end

    step_A!(r, rtmp, dt/2, v)

    return nothing
end

function propagate_mapping_variables!(cache, calc, X, P, rtmp, dt)

    @unpack C, D, tmp_vec1, tmp_vec2 = cache

    eigen = Calculators.get_eigen(calc, rtmp)

    set_C_propagator!(C, cache, eigen, dt)
    set_D_propagator!(D, cache, eigen, dt)

    # tmp_vec1 = C*X - D*P
    mul!(tmp_vec1, C, X)
    mul!(tmp_vec1, D, P, -1.0, 1.0)

    # tmp_vec2 = C*P + D*X
    mul!(tmp_vec2, C, P)
    mul!(tmp_vec2, D, X, 1.0, 1.0)

    copyto!(X, tmp_vec1)
    copyto!(P, tmp_vec2)
end

"Get the `C` propagator for the mapping variables."
function set_C_propagator!(C, cache, eigen::Eigen, dt::Real)
    fill!(C, zero(eltype(C)))
    for i in axes(C,1)
        C[i,i] = cos(eigen.values[i] * dt)
    end
    transform_matrix!(C, eigen.vectors, cache.tmp_mat)
end

"Get the `D` propagator for the mapping variables."
function set_D_propagator!(D, cache, eigen::Eigen, dt::Real)
    fill!(D, zero(eltype(D)))
    for i in axes(D,1)
        D[i,i] = sin(-eigen.values[i] * dt)
    end
    transform_matrix!(D, eigen.vectors, cache.tmp_mat)
end

function transform_matrix!(M, transform, tmp_mat)
    mul!(tmp_mat, M, transform')
    mul!(M, transform, tmp_mat)
    return nothing
end

"Get the `Γ` variable used to calculate the nuclear propagators."
function set_gamma!(Γ::AbstractMatrix, W::AbstractMatrix, Λ::AbstractVector, dt::Real)
    for i in axes(Γ, 2)
        Γ[i,i] = W[i,i] * dt
        for j in axes(Γ, 1)
            if j != i
                Γ[j,i] = sin((Λ[i] - Λ[j])*dt) * W[j,i] / (Λ[i] - Λ[j])
            end
        end
    end
    return Γ
end

"Get the `Ξ` variable used to calculate the nuclear propagators."
function set_xi!(Ξ::AbstractMatrix, W::AbstractMatrix, Λ::AbstractVector, dt::Real)
    for i in axes(Ξ, 2)
        for j in axes(Ξ, 1)
            if j != i
                Ξ[j,i] = (1 - cos((Λ[i] - Λ[j])*dt)) * W[j,i] / (Λ[i] - Λ[j])
            end
        end
    end
    return Ξ
end

"""
Get the force due to the mapping variables.

Equivalent to this but doesn't allocate: 
    return 0.5 * (q'*E*q + p'*E*p) - q'*F*p
"""
function get_mapping_nuclear_force(q::AbstractVector, p::AbstractVector,
                           E::AbstractMatrix, F::AbstractMatrix, tmp_vec)
    force = zero(eltype(q))
    mul!(tmp_vec, E, q)
    force += dot(q, tmp_vec)
    mul!(tmp_vec, E, p)
    force += dot(p, tmp_vec)
    mul!(tmp_vec, F, p)
    force -= 2dot(q, tmp_vec)
    return force / 2
end

DynamicsMethods.select_algorithm(::Simulation{<:MappingVariableMethods.SpinMappingW}) = MInt()
