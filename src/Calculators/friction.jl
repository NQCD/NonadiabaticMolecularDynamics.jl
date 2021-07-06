
mutable struct FrictionCalculator{T,M} <: AbstractFrictionCalculator{M}
    model::M
    potential::T
    derivative::Matrix{T}
    friction::Matrix{T}
    function FrictionCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        new{T,M}(model, 0, zeros(DoFs, atoms), zeros(DoFs*atoms, DoFs*atoms))
    end
end

struct RingPolymerFrictionCalculator{T,M} <: AbstractFrictionCalculator{M}
    model::M
    potential::Vector{T}
    derivative::Array{T,3}
    friction::Array{T,3}
    function RingPolymerFrictionCalculator{T}(model::M, DoFs::Integer, atoms::Integer, beads::Integer) where {T,M<:Model}
        new{T,M}(model, zeros(beads), zeros(DoFs, atoms, beads), zeros(DoFs*atoms, DoFs*atoms, beads))
    end
end

struct DiabaticFrictionCalculator{T,M} <: AbstractDiabaticCalculator{M}
    model::M
    potential::Hermitian{T,Matrix{T}}
    derivative::Matrix{Hermitian{T,Matrix{T}}}
    eigenvalues::Vector{T}
    eigenvectors::Matrix{T}
    adiabatic_derivative::Matrix{Matrix{T}}
    nonadiabatic_coupling::Matrix{Matrix{T}}
    friction::Matrix{T}
    tmp_mat::Matrix{T}
    tmp_mat_complex1::Matrix{Complex{T}}
    tmp_mat_complex2::Matrix{Complex{T}}
    function DiabaticFrictionCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        potential = Hermitian(zeros(model.n_states, model.n_states))
        derivative = [Hermitian(zeros(model.n_states, model.n_states)) for i=1:DoFs, j=1:atoms]
        eigenvalues = zeros(model.n_states)
        eigenvectors = zeros(model.n_states, model.n_states)
        adiabatic_derivative = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms]
        nonadiabatic_coupling = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms]
        friction = zeros(DoFs*atoms, DoFs*atoms)
        tmp_mat = zeros(T, model.n_states, model.n_states)
        tmp_mat_complex1 = zeros(Complex{T}, model.n_states, model.n_states)
        tmp_mat_complex2 = zeros(Complex{T}, model.n_states, model.n_states)
        new{T,M}(model, potential, derivative, eigenvalues, eigenvectors,
                 adiabatic_derivative, nonadiabatic_coupling, friction,
                 tmp_mat, tmp_mat_complex1, tmp_mat_complex2)
    end
end
function Calculator(model::AdiabaticFrictionModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    FrictionCalculator{T}(model, DoFs, atoms)
end
function Calculator(model::DiabaticFrictionModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    DiabaticFrictionCalculator{T}(model, DoFs, atoms)
end
function Calculator(model::AdiabaticFrictionModel, DoFs::Integer, atoms::Integer, beads::Integer, T::Type=Float64)
    RingPolymerFrictionCalculator{T}(model, DoFs, atoms, beads)
end

function evaluate_friction!(calc::AbstractFrictionCalculator, R::AbstractMatrix)
    friction!(calc.model, calc.friction, R)
end

function evaluate_friction!(calc::AbstractFrictionCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        friction!(calc.model, calc.friction[:,:,i], R[:,:,i])
    end
end

@doc raw"""
    evaluate_friction!(calc::DiabaticFrictionCalculator, R::AbstractMatrix)

Evaluate the electronic friction for a model given in the diabatic representation.

Requires that `adiabatic_derivative` and `eigenvalues` be precomputed.

```math
γ = 2πħ ∑ⱼ <1|dH|j><j|dH|1> δ(ωⱼ) / ωⱼ
```
Note that the delta function is approximated by a normalised gaussian.
"""
function evaluate_friction!(calc::DiabaticFrictionCalculator, R::AbstractMatrix)

    gauss(x, σ) = exp(-0.5 * x^2 / σ^2) / (σ*sqrt(2π))

    DoFs = size(R)[1]
    calc.friction .= 0
    for i in axes(R, 2) # Atoms
        for j in axes(R, 1) # DoFs
            for m=2:calc.model.n_states
                ω = calc.eigenvalues[m] - calc.eigenvalues[1]
                g = gauss(ω, calc.model.σ) / ω
                calc.friction[j+(i-1)*DoFs] += 2π*abs2(calc.adiabatic_derivative[j,i][m,1])*g
            end
        end
    end
end