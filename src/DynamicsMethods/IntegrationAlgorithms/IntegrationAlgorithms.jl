
module IntegrationAlgorithms

using NQCDynamics:
    NQCDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsMethods,
    DynamicsUtils,
    Calculators,
    natoms, nbeads, ndofs

using OrdinaryDiffEq: OrdinaryDiffEq, OrdinaryDiffEqAlgorithm
using StochasticDiffEq: StochasticDiffEq

struct BCB <: OrdinaryDiffEqAlgorithm end
struct BCBFull <: OrdinaryDiffEqAlgorithm end

struct BCBwithTsit5{T<:OrdinaryDiffEqAlgorithm} <: OrdinaryDiffEqAlgorithm
    electronic_algorithm::T
end

struct BCBWavefunction <: OrdinaryDiffEqAlgorithm end

struct BABwithTsit5{T<:OrdinaryDiffEqAlgorithm} <: OrdinaryDiffEqAlgorithm
    electronic_algorithm::T
end

"""
    RingPolymerMInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm

Second order symplectic momentum integral algorithm applied to NRPMD.

# Reference

[J. Chem. Phys. 148, 102326 (2018)](https://doi.org/10.1063/1.5005557)
"""
struct RingPolymerMInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

"""
    MInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm

Second order symplectic momentum integral algorithm.

# Reference

[J. Chem. Phys. 148, 102326 (2018)](https://doi.org/10.1063/1.5005557)
"""
struct MInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end
struct VerletwithElectronics <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

struct MDEF_BAOAB <: StochasticDiffEq.StochasticDiffEqAlgorithm end
struct BCOCB <: StochasticDiffEq.StochasticDiffEqAlgorithm end


"""
This block of code below is to select the algorithm for each type of simulation.

    

"""

DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping}) = BCBwithTsit5(OrdinaryDiffEq.Tsit5())
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.EhrenfestMethods.AbstractEhrenfest}) = BCBwithTsit5(OrdinaryDiffEq.Tsit5())
DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.ClassicalMethods.AbstractMDEF}) = MDEF_BAOAB()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.AbstractMDEF}) = BCOCB()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.ThermalLangevin}) = BCOCB()
DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.MappingVariableMethods.SpinMappingW}) = MInt()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.MappingVariableMethods.NRPMD}) = RingPolymerMInt()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.MappingVariableMethods.eCMM}) = RingPolymerMInt()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{DynamicsMethods.ClassicalMethods.Classical}) = BCB()
DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.SurfaceHoppingMethods.AbstractIESH}) = VerletwithElectronics()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.SurfaceHoppingMethods.AbstractIESH}) = BCBWavefunction()
DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.EhrenfestMethods.Ehrenfest}) = BABwithTsit5(OrdinaryDiffEq.Tsit5())
DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.EhrenfestMethods.EhrenfestNA}) = VerletwithElectronics()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.EhrenfestMethods.EhrenfestNA}) = BCBWavefunction()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.SurfaceHoppingMethods.ClassicalMasterEquation}) = BCBFull()

export BCB
export BCBwithTsit5
export BABwithTsit5
export RingPolymerMInt
export MInt
export MDEF_BAOAB
export BCOCB
export BCBWavefunction
export BCBFull

include("mdef_baoab.jl")
include("bcocb.jl")
include("mint.jl")
include("ringpolymer_mint.jl")
include("bcb_electronics.jl")
include("bab_electronics.jl")
include("bcb.jl")
include("steps.jl")
include("verlet_with_electronics.jl")
include("bcb_wavefunction.jl")
include("bcb_full.jl")

end # module
