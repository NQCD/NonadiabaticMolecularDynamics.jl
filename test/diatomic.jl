using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions

atoms = Atoms([:N, :O])
model = Models.DiatomicHarmonic()
sim = Simulation(atoms, model, Dynamics.Classical())

@testset "calculate_diatomic_energy" begin
    @test Diatomic.calculate_diatomic_energy(sim, 1.0, normal_vector=rand(3)) ≈ 0.0 atol=1e-3
    @test Diatomic.calculate_diatomic_energy(sim, 2.0, normal_vector=rand(3)) ≈ 0.5
end

@testset "calculate_force_constant" begin
    @test Diatomic.calculate_force_constant(sim, 1.0) ≈ 1
    model2 = Models.DiatomicHarmonic(r₀=6.3)
    sim2 = Simulation(atoms, model2, Dynamics.Classical())
    @test Diatomic.calculate_force_constant(sim2, 6.3) ≈ 1
end
