using NQCModels: AdiabaticModels.AdiabaticChemicalEnvironmentMLIP

function evaluate_potential!(calc::RingPolymerAdiabaticCalculator{T, M}, R::AbstractArray{S, 3}) where {T, M::AdiabaticModels.AdiabaticChemicalEnvironmentMLIP, S}
    calc.stats[:potential] += 1
    potential_vector = @views [R[:, :, i] for i in axes(R, 3)]
    NQCModels.potential(calc.model, calc.model.atoms, potential_vector, calc.model.cell)
    return nothing
end

function evaluate_derivative!(calc::RingPolymerAdiabaticCalculator{T, M}, R::AbstractArray{S, 3}) where {T, M::AdiabaticModels.AdiabaticChemicalEnvironmentMLIP, S}
    calc.stats[:derivative] += 1
    derivative_vector = @views [R[:, :, i] for i in axes(R, 3)]
    NQCModels.derivative(calc.model, calc.model.atoms, derivative_vector, calc.model.cell)
    return nothing
end

function evaluate_potential!(calc::RingPolymerFrictionCalculator{T, M}, R::AbstractArray{S, 3}) where {T, M::AdiabaticModels.AdiabaticChemicalEnvironmentMLIP, S}
    calc.stats[:potential] += 1
    potential_vector = @views [R[:, :, i] for i in axes(R, 3)]
    NQCModels.potential(calc.model, calc.model.atoms, potential_vector, calc.model.cell)
    return nothing
end

function evaluate_derivative!(calc::RingPolymerFrictionCalculator{T, M}, R::AbstractArray{S, 3}) where {T, M::AdiabaticModels.AdiabaticChemicalEnvironmentMLIP, S}
    calc.stats[:derivative] += 1
    derivative_vector = @views [R[:, :, i] for i in axes(R, 3)]
    NQCModels.derivative(calc.model, calc.model.atoms, derivative_vector, calc.model.cell)
    return nothing
end
