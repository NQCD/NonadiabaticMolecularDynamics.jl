"""
    Calculators

This module exists to bridge the gap between the `Models` and the `Dynamics`.

Here we provide functions and types for evaluating and storing quantities obtained from the
`Models`.
In addition any further manipulation of those quantities, such as computing eigenvalues,
is included here.

This module is largely needed to facilitate seemless integration of both ring polymer and
classical dynamics, using the same models and functions for both.
Specific ring polymer types are provided that have the extra fields and methods needed
to evaluate the quantities for each bead. 
"""
module Calculators

using LinearAlgebra: LinearAlgebra, Hermitian, I
using StaticArrays: SMatrix, SVector

using NonadiabaticModels: NonadiabaticModels, Model, nstates
using NonadiabaticModels.AdiabaticModels: AdiabaticModel
using NonadiabaticModels.DiabaticModels: DiabaticModel, DiabaticFrictionModel
using NonadiabaticModels.FrictionModels: AdiabaticFrictionModel

using NonadiabaticMolecularDynamics: RingPolymers

"""
    AbstractCalculator{M<:Model}

Top-level type for all calculators.

Each concrete calculator contains the `Model` and the fields to store the quantities
obtained from the model.
"""
abstract type AbstractCalculator{M<:Model} end
abstract type AbstractAdiabaticCalculator{M<:AdiabaticModel} <: AbstractCalculator{M} end
abstract type AbstractDiabaticCalculator{M<:Union{DiabaticFrictionModel,DiabaticModel}} <: AbstractCalculator{M} end
abstract type AbstractStaticDiabaticCalculator{M} <: AbstractDiabaticCalculator{M} end
abstract type AbstractFrictionCalculator{M<:AdiabaticFrictionModel} <: AbstractCalculator{M} end

mutable struct AdiabaticCalculator{T,M} <: AbstractAdiabaticCalculator{M}
    model::M
    potential::T
    derivative::Matrix{T}
    function AdiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        new{T,M}(model, 0, zeros(DoFs, atoms))
    end
end

struct RingPolymerAdiabaticCalculator{T,M} <: AbstractAdiabaticCalculator{M}
    model::M
    potential::Vector{T}
    derivative::Array{T,3}
    function RingPolymerAdiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer, beads::Integer) where {T,M<:Model}
        new{T,M}(model, zeros(beads), zeros(DoFs, atoms, beads))
    end
end

mutable struct DiabaticCalculator{T,M,S,L} <: AbstractDiabaticCalculator{M}
    model::M
    potential::Hermitian{T,SMatrix{S,S,T,L}}
    derivative::Matrix{Hermitian{T,SMatrix{S,S,T,L}}}
    eigenvalues::SVector{S,T}
    eigenvectors::SMatrix{S,S,T,L}
    adiabatic_derivative::Matrix{SMatrix{S,S,T,L}}
    nonadiabatic_coupling::Matrix{SMatrix{S,S,T,L}}
    tmp_mat::Matrix{T}
    tmp_mat_complex1::Matrix{Complex{T}}
    tmp_mat_complex2::Matrix{Complex{T}}
    function DiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        n = nstates(model)
        matrix_template = NonadiabaticModels.DiabaticModels.matrix_template(model, T)
        vector_template = NonadiabaticModels.DiabaticModels.vector_template(model, T)

        potential = Hermitian(matrix_template)
        derivative = [Hermitian(matrix_template) for _=1:DoFs, _=1:atoms]
        eigenvalues = vector_template
        eigenvectors = matrix_template + I
        adiabatic_derivative = [matrix_template for _ in CartesianIndices(derivative)]
        nonadiabatic_coupling = [matrix_template for _ in CartesianIndices(derivative)]
        tmp_mat = zeros(T, n, n)
        tmp_mat_complex1 = zeros(Complex{T}, n, n)
        tmp_mat_complex2 = zeros(Complex{T}, n, n)

        new{T,M,n,n^2}(model,
            potential, derivative, eigenvalues, eigenvectors,
            adiabatic_derivative, nonadiabatic_coupling,
            tmp_mat, tmp_mat_complex1, tmp_mat_complex2)
    end
end

mutable struct RingPolymerDiabaticCalculator{T,M,S,L} <: AbstractDiabaticCalculator{M}
    model::M
    potential::Vector{Hermitian{T,SMatrix{S,S,T,L}}}
    derivative::Array{Hermitian{T,SMatrix{S,S,T,L}},3}
    eigenvalues::Vector{SVector{S,T}}
    eigenvectors::Vector{SMatrix{S,S,T,L}}
    adiabatic_derivative::Array{SMatrix{S,S,T,L},3}
    nonadiabatic_coupling::Array{SMatrix{S,S,T,L},3}

    centroid_potential::Hermitian{T,SMatrix{S,S,T,L}}
    centroid_derivative::Matrix{Hermitian{T,SMatrix{S,S,T,L}}}
    centroid_eigenvalues::SVector{S,T}
    centroid_eigenvectors::SMatrix{S,S,T,L}
    centroid_adiabatic_derivative::Matrix{SMatrix{S,S,T,L}}
    centroid_nonadiabatic_coupling::Matrix{SMatrix{S,S,T,L}}

    tmp_mat::Matrix{T}
    tmp_mat_complex1::Matrix{Complex{T}}
    tmp_mat_complex2::Matrix{Complex{T}}
    function RingPolymerDiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer, beads::Integer) where {T,M<:Model}
        n = nstates(model)
        matrix_template = NonadiabaticModels.DiabaticModels.matrix_template(model, T)
        vector_template = NonadiabaticModels.DiabaticModels.vector_template(model, T)

        potential = [Hermitian(matrix_template) for _=1:beads]
        derivative = [Hermitian(matrix_template) for _=1:DoFs, _=1:atoms, _=1:beads]
        eigenvalues = [vector_template for _=1:beads]
        eigenvectors = [matrix_template + I for _=1:beads]
        adiabatic_derivative = [matrix_template for _=1:DoFs, _=1:atoms, _=1:beads]
        nonadiabatic_coupling = [matrix_template for _=1:DoFs, _=1:atoms, _=1:beads]

        centroid_potential = Hermitian(matrix_template)
        centroid_derivative = [Hermitian(matrix_template) for _=1:DoFs, _=1:atoms]
        centroid_eigenvalues = vector_template
        centroid_eigenvectors = matrix_template + I
        centroid_adiabatic_derivative = [matrix_template for _ in CartesianIndices(centroid_derivative)]
        centroid_nonadiabatic_coupling = [matrix_template for _ in CartesianIndices(centroid_derivative)]

        tmp_mat = zeros(T, n, n)
        tmp_mat_complex1 = zeros(Complex{T}, n, n)
        tmp_mat_complex2 = zeros(Complex{T}, n, n)
        new{T,M,n,n^2}(model, potential, derivative, eigenvalues, eigenvectors, adiabatic_derivative, nonadiabatic_coupling,
            centroid_potential, centroid_derivative, centroid_eigenvalues, centroid_eigenvectors, 
            centroid_adiabatic_derivative, centroid_nonadiabatic_coupling,
            tmp_mat, tmp_mat_complex1, tmp_mat_complex2)
    end
end

function Calculator(model::DiabaticModel, DoFs::Integer, atoms::Integer, t::Type{T}) where {T}
    DiabaticCalculator{t}(model, DoFs, atoms)
end
function Calculator(model::AdiabaticModel, DoFs::Integer, atoms::Integer, t::Type{T}) where {T}
    AdiabaticCalculator{t}(model, DoFs, atoms)
end
function Calculator(model::DiabaticModel, DoFs::Integer, atoms::Integer, beads::Integer, t::Type{T}) where {T}
    RingPolymerDiabaticCalculator{t}(model, DoFs, atoms, beads)
end
function Calculator(model::AdiabaticModel, DoFs::Integer, atoms::Integer, beads::Integer, t::Type{T}) where {T}
    RingPolymerAdiabaticCalculator{t}(model, DoFs, atoms, beads)
end

function evaluate_potential!(calc::AbstractCalculator, R)
    calc.potential = NonadiabaticModels.potential(calc.model, R)
end

function evaluate_potential!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        calc.potential[i] = NonadiabaticModels.potential(calc.model, R[:,:,i])
    end
end

function evaluate_centroid_potential!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    calc.centroid_potential = NonadiabaticModels.potential(calc.model, RingPolymers.get_centroid(R))
end

function evaluate_derivative!(calc::AbstractCalculator, R)
    NonadiabaticModels.derivative!(calc.model, calc.derivative, R)
end

function evaluate_derivative!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        NonadiabaticModels.derivative!(calc.model, calc.derivative[:,:,i], R[:,:,i])
    end
end

function evaluate_centroid_derivative!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    NonadiabaticModels.derivative!(calc.model, calc.centroid_derivative, RingPolymers.get_centroid(R))
end

function eigen!(calc::DiabaticCalculator)
    eig = LinearAlgebra.eigen(calc.potential)
    calc.eigenvalues = eig.values
    calc.eigenvectors = correct_phase(eig.vectors, calc.eigenvectors)
    return nothing
end

function centroid_eigen!(calc::RingPolymerDiabaticCalculator)
    eig = LinearAlgebra.eigen(calc.centroid_potential)
    calc.centroid_eigenvalues = eig.values
    calc.centroid_eigenvectors = correct_phase(eig.vectors, calc.centroid_eigenvectors)
    return nothing
end

function eigen!(calc::RingPolymerDiabaticCalculator)
    for i=1:length(calc.potential)
        eig = LinearAlgebra.eigen(calc.potential[i])
        calc.eigenvalues[i] = eig.values
        calc.eigenvectors[i] = correct_phase(eig.vectors, calc.eigenvectors[i])
    end
    return nothing
end

function correct_phase(new_vectors::SMatrix, old_vectors::SMatrix)
    n = size(new_vectors, 1)
    vect = SVector{n}(sign(LinearAlgebra.dot(new_vectors[:,i], old_vectors[:,i])) for i=1:n)
    return new_vectors .* vect'
end

function transform_derivative!(calc::AbstractDiabaticCalculator)
    for I in eachindex(calc.derivative)
        calc.adiabatic_derivative[I] = calc.eigenvectors' * calc.derivative[I] * calc.eigenvectors
    end
end

function transform_centroid_derivative!(calc::RingPolymerDiabaticCalculator)
    for I in eachindex(calc.centroid_derivative)
        calc.centroid_adiabatic_derivative[I] = calc.centroid_eigenvectors' * calc.centroid_derivative[I] * calc.centroid_eigenvectors
    end
end

function transform_derivative!(calc::RingPolymerDiabaticCalculator)
    for i in axes(calc.derivative, 3) # Beads
        for j in axes(calc.derivative, 2) # Atoms
            for k in axes(calc.derivative, 1) # DoFs
                calc.adiabatic_derivative[k,j,i] = calc.eigenvectors[i]' * calc.derivative[k,j,i] * calc.eigenvectors[i]
            end
        end
    end
end

function evaluate_nonadiabatic_coupling!(calc::AbstractDiabaticCalculator)
    for I in eachindex(calc.adiabatic_derivative)
        calc.nonadiabatic_coupling[I] = evaluate_nonadiabatic_coupling(calc.adiabatic_derivative[I], calc.eigenvalues)
    end
end

function evaluate_centroid_nonadiabatic_coupling!(calc::RingPolymerDiabaticCalculator)
    for I in eachindex(calc.centroid_adiabatic_derivative)
        calc.centroid_nonadiabatic_coupling[I] = evaluate_nonadiabatic_coupling(calc.centroid_adiabatic_derivative[I], calc.centroid_eigenvalues)
    end
end

function evaluate_nonadiabatic_coupling!(calc::RingPolymerDiabaticCalculator)
    for i in axes(calc.nonadiabatic_coupling, 3) # Beads
        for I in CartesianIndices(size(calc.adiabatic_derivative)[1:2])
            calc.nonadiabatic_coupling[I,i] = evaluate_nonadiabatic_coupling(calc.adiabatic_derivative[I,i], calc.eigenvalues[i])
        end
    end
end

"""
# References

- HammesSchifferTully_JChemPhys_101_4657_1994 Eq. (32)
- SubotnikBellonzi_AnnuRevPhyschem_67_387_2016, section 2.3
"""
function evaluate_nonadiabatic_coupling(adiabatic_derivative::SMatrix, eigenvalues::SVector)
    n = length(eigenvalues)
    SMatrix{n,n}(
        (i != j ? adiabatic_derivative[j,i] / (eigenvalues[i] - eigenvalues[j]) : 0
        for j=1:n, i=1:n))
end

"""
Evaluates all electronic properties for the current position `r`.

# Properties evaluated:
- Diabatic potential
- Diabatic derivative
- Eigenvalues and eigenvectors
- Adiabatic derivative
- Nonadiabatic coupling
"""
function update_electronics!(calculator::AbstractDiabaticCalculator, r::AbstractArray)
    evaluate_potential!(calculator, r)
    evaluate_derivative!(calculator, r)
    eigen!(calculator)
    transform_derivative!(calculator)
    evaluate_nonadiabatic_coupling!(calculator)
end

function update_electronics!(calculator::RingPolymerDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    evaluate_potential!(calculator, r)
    evaluate_derivative!(calculator, r)
    eigen!(calculator)
    transform_derivative!(calculator)
    evaluate_nonadiabatic_coupling!(calculator)

    update_centroid_electronics!(calculator, r)
end

function update_centroid_electronics!(calculator::RingPolymerDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    evaluate_centroid_potential!(calculator, r)
    evaluate_centroid_derivative!(calculator, r)
    centroid_eigen!(calculator)
    transform_centroid_derivative!(calculator)
    evaluate_centroid_nonadiabatic_coupling!(calculator)
end

include("large_diabatic.jl")
include("friction.jl")

end # module
