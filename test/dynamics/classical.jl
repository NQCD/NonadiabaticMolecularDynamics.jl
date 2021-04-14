using NonadiabaticMolecularDynamics
using Test
using Unitful

include("utils.jl")

atoms = Atoms([:H])
model = Models.Harmonic()
sim = Simulation{Classical}(atoms, model; DoFs=1)

v = rand(1, length(atoms))
r = rand(1, length(atoms))
dr = zero(r)
dv = zero(v)

test_velocity!(sim)
test_acceleration!(sim)
test_motion!(sim)

u0 = ClassicalDynamicals(v, r)
sol = Dynamics.run_trajectory(u0, (0.0, 1000.0), sim; dt=0.1)

e0 = evaluate_hamiltonian(sim, ClassicalDynamicals(sol.u[1]))
e1 = evaluate_hamiltonian(sim, ClassicalDynamicals(sol.u[end]))
@test e0 ≈ e1 rtol=1e-2

sim = RingPolymerSimulation{Classical}(atoms, model, 10; DoFs=1, temperature=100u"K")

v = RingPolymerArray(rand(1, length(atoms), length(sim.beads)))
r = RingPolymerArray(rand(1, length(atoms), length(sim.beads)))
dv = zero(v)
dr = zero(r)

u0 = ClassicalDynamicals(v, r)
sol = Dynamics.run_trajectory(u0, (0.0, 1000.0), sim; dt=0.1)

e0 = evaluate_hamiltonian(sim, ClassicalDynamicals(sol.u[1]))
e1 = evaluate_hamiltonian(sim, ClassicalDynamicals(sol.u[end]))
@test e0 ≈ e1 rtol=1e-2