```@setup logging
@info "Expanding src/examples/reactive_scattering.md..."
start_time = time()
```
# [Reactive scattering from a metal surface](@id example-h2scattering)

Our implementation allows us to simulate vibrational de-excitation probability during reactive scattering events at metal surfaces for any diatomic molecule 
with a suitable model to describe energies and forces (and friction coefficients for MDEF simulations). 
Here, we investigate the reactive scattering of hydrogen on a Cu(111) metal surface as a prototypical example.

To run this kind of simulation, a set of initial positions and velocities (``\mathbf{R}`` and ``\mathbf{\dot{R}}``) with
ro-vibrational quantum states ``\nu`` and ``j`` have to be generated
(see [EBK quantisation](@ref ebk-sampling)).
With a specific ro-vibrational quantum state it is possible to compute different properties
after molecular collision and energy transfer with the metal surface like the vibrational
de-excitation probabilities discussed here.

In order to reproduce the state-to-state vibrational de-excitation probability results presented originally by [Maurer2019](@cite) for this system, the same initial conditions were generated with [`QuantisedDiatomic.generate_configurations`](@ref InitialConditions.QuantisedDiatomic.generate_configurations)
setting the initial ro-vibrational quantum state to (``\nu=2, j=0``) as was explored in the original paper. 

As shown earlier in the [EBK documentation](@ref ebk-sampling) we are able to generate
a semiclassically quantised distribution for a diatomic molecule on a collision course
with a metal surface.
In this example we follow the [EBK example](@ref ebk-sampling) using [`machine learning potential`](@ref ml-pes-models)
to prepare our initial distribution and run our simulation.

Specifically, we have produced a set of initial conditions with different translational energy (`translational_energy` keyword) ranging from 0.2 to 1.4 eV, locating the hydrogen molecule 8 Å away from the metal
surface (`height` keyword) with a normal incidence.

!!! note "Atomic units"

    As usual, all quantities default to atomic units. Here we use [Unitful](https://painterqubits.github.io/Unitful.jl/stable/)
    to input the translational energy and height using different units, where they are later converted internally.

```@example h2scatter
using NQCDynamics
using Unitful
using NQCDynamics.InitialConditions: QuantisedDiatomic
using JLD2
using PyCall

mace_calc = pyimport("mace.calculators")

atoms = Atoms([:H, :H])
cell = PeriodicCell([11.1175 -5.5588 0.0; 0.0 9.628 0.0; 0.0 0.0 70.3079])
positions = [0.0 0.0; 0.0 0.0; 0.0 0.73]
atoms_ase = NQCDynamics.convert_from_ase_atoms(atoms, positions, cell)
ids_adsorbate = [55,56]

# load PES model
calculator = mace_calc.MACECalculator(model_path="../assets/mace/h2cu.model", device="cpu", default_dtype="float32")
atoms_ase.set_calculator(calculator)
model = AdiabaticASEModel(atoms_ase)

sim = Simulation(atoms, model; cell=cell)

ν, J = 2, 0     # selected ro-vibrational quantum states  
nsamples = 10  # number of configurations      
Ek = 0.5u"eV"   # Translational energy [eV] ; range considered [0.2-1.4] eV
z = 8.0u"Å"     # Height [Å]  ; fixed at 8 Å
z = ang_to_au(z)

configurations = QuantisedDiatomic.generate_configurations(sim, ν, J;
    samples=nsamples, translational_energy=Ek, height=z)
v_h2 = first.(configurations)
r_h2 = last.(configurations)


# Re build the position and velocities for the whole system Cu56H2
dof = 3
ase_io = pyimport("ase.io")
surface_ase = ase_io.read("../assets/h2cu/h2cu_surf.traj@:nsamples")

n_atoms = length(atoms)
v = [zeros(dof,n_atoms) for i=1:length(r_h2)]
r = [zeros(dof,n_atoms) for i=1:length(r_h2)]

for i in 1:length(r_h2)
    surf_atoms, surf_positions, surf_cell = NQCDynamics.convert_from_ase_atoms(surface_ase[i])
    r[i][:,1:n_atoms-2] .= surf_positions
    r[i][:,n_atoms-1:n_atoms] .= r_h2[i]  

    v[i][:,1:n_atoms-2] .= austrip.(transpose(surface_ase[i].get_velocities().*ase_units.fs)*u"Å/fs")
    v[i][:,n_atoms-1:n_atoms] .= v_h2[i]
end

atoms_ase = ase_io.read("../assets/h2cu/h2cu_init.in")
atoms, positions, cell = convert_from_ase_atoms(atoms_ase)

distribution = DynamicalDistribution(v, r, (dof,length(atoms_ase)))
nothing # hide
```

!!! tip "Saving the distribution"

    Generally it will be desirable to generate a distribution once and re-use it for multiple dynamics simulations.
    The simplest way to do this is to save the distribution using [JLD2.jl](https://juliaio.github.io/JLD2.jl/dev/).
    Refer to [Saving and loading](@ref saving-and-loading) to learn more.

In order to produce an unweighted distribution, the lateral and angular orientation are randomly selected within the unit cell.
As an example of the spacial and orientation distribution generated with this module, a subset of data (300 configurations) is shown below.
To run our production simulations, however, a set of 80,000 initial velocities and positions were used.

![initial conditions](../assets/figures/icond_scatter.png)

## Data analysis and truncation function

Since we are interested in the dynamics only when the molecule is close to the surface,
we can use a callback to terminate the simulation early to save us some time.
This requires defining a function that returns `true` when we want the simulation to
terminate.
This means we can set our time span relatively long since we expect most simulations to
terminate before reaching the time limit.

```@example h2scatter
using Statistics: mean
using LinearAlgebra: norm

h2distance(p) = norm(p[:,ids_adsorbate[1]] .- p[:,ids_adsorbate[2]])

"Terminates simulation if returns `true`." 
mutable struct TrajectoryTerminator
    h2_indices
    ads_height_cutoff
    ads_dist_cutoff
    n_atoms_layer
end
function (t::TrajectoryTerminator)(u, t, integrator)::Bool
    R = get_positions(u)
    com_h2_z = minimum(R[3,t.h2_indices[1]:t.h2_indices[2]])
    top_surface_avg_z = mean(R[3,end-Int(t.n_atoms_layer)-1:end-2])
    zcom = au_to_ang(com_h2_z-top_surface_avg_z) # Convert vertical centre of mass to angstrom
    if zcom > t.ads_height_cutoff                         # Scattering event
        return true
    elseif au_to_ang(h2distance(R)) > t.ads_dist_cutoff   # Reactive event
        return true
    else
        return false
    end
end

termination_condition = TrajectoryTerminator(ids_adsorbate, 8.1, 2.5, 9)
terminate = DynamicsUtils.TerminatingCallback(termination_condition)
tspan = (0.0, 420.0u"fs")
nothing # hide
```
In this example, we consider the outcome a reactive event if the H-H
bond length is larger than 2.5 Å in any point of during the trajectory and a
scattering event if the molecule rebounds to a vertical distance from the metal
surface greater than 8.1 Å.

## MDEF with the LDFA

Now that we have set up the initial distribution and some of our simulation parameters,
we can choose which form of friction we would like use.
First, let's use the cube-based density implementation for LDFA friction provided by the
[FrictionProviders.jl](@ref friction-providers).
This takes a `.cube` file containing the electron density and will provide the friction
during the dynamics.
Here we initialize the MACE-based interatomic potential, together with the cube-based calculator.

```@example h2scatter
# load PES model
calculator = mace_calc.MACECalculator(model_path="../assets/mace/h2cu.model", device="cpu", default_dtype="float32")
atoms_ase.set_calculator(calculator)
model_pes = AdiabaticASEModel(atoms_ase)

# load cube EFT calculator
using FrictionProviders
density_model = CubeLDFA("../assets/friction/test.cube", cell)
model_eft = LDFAFriction(density_model, atoms; friction_atoms=ids_adsorbate)
```

Now we can pass all the variables defined so far to the `Simulation` and run multiple
trajectories using [`run_dynamics`](@ref).

```julia
model = CompositeFrictionModel(model_pes, model_eft)
sim = Simulation{MDEF}(atoms, model, cell=cell, temperature=300u"K")
ensemble = run_dynamics(sim, tspan, distribution; selection=1:nsamples,
    dt=0.1u"fs", output=OutputPosition, trajectories=nsamples, callback=terminate)
```



## MDEF with machine-learned LDFA 

Above, we used the cube-based LDFA interpretation of MDEF to perform the simulation.
However, in this scheme, surface atoms have to stay fixed to match the cube densities. To include surface temperature effects machine learning models can be trained that predict densities at any surface configuration.
Here, we run MDEF simulation employing two machine learning models, to predict adiabatic PES (based on [MACE](https://github.com/ACEsuit/mace)) and surface electron density (based on [ACEpotentials.jl](https://github.com/ACEsuit/ACEpotentials.jl)), utilizing [FrictionProviders.jl](https://github.com/NQCD/FrictionProviders.jl), to run MDEF simulation.

```julia
using ACEpotentials
ace_model, ace_model_meta = ACEpotentials.load_model("../assets/ace_ldfa/model.json")
density_model = AdiabaticModels.ACEpotentialsModel(atoms, cell, ace_model) 
ace_density_model = AceLDFA(density_model)
model_eft = LDFAFriction(ace_density_model, atoms; friction_atoms=ids_adsorbate)

model = CompositeFrictionModel(model_pes, model_eft)

sim = Simulation{MDEF}(atoms, model, cell=cell, temperature=300u"K")
ensemble = run_dynamics(sim, tspan, distribution; selection=1:nsamples,
    dt=0.1u"fs", output=OutputPosition, trajectories=nsamples, callback=terminate)
```

## MDEF with machine-learned TDPT (ODF)

Above, we used the LDFA interpretation of MDEF to perform the simulation. However, an alternative, TDPT (otherwise known as ODF) method can be used that provides full friction tensor. TDPT ML models can be incorporated in our simulation in a similar way as LDFA models, through (FrictionProviders.jl)[https://github.com/NQCD/FrictionProviders.jl]. Here, we show how this can be done for (ACEds)[https://github.com/ACEsuit/ACEds.jl] models. Such models require the usage of (JuLIP)[https://github.com/JuliaMolSim/JuLIP.jl]-type atoms. We can easily convert our ASE-type atoms into such format:

```julia
using ASE, JuLIP
atoms_ase_jl = ASE.ASEAtoms(atoms_ase)
atoms_julip = JuLIP.Atoms(atoms_ase_jl)
```

Having our atoms, ACEds EFT models can be then initialized and combined with previously loaded potential model.

```julia
using ACE
using ACEds.FrictionModels
using ACEds.FrictionModels: Gamma
aceds_model = ACEdsODF(read_dict(load_dict("../assets/ace_friction/eft_ac.model")), Gamma, atoms_julip)
model_eft = ODFriction(ace_model; friction_atoms=ids_adsorbate)

model = CompositeFrictionModel(model_pes, model_eft)

sim = Simulation{MDEF}(atoms, model, cell=cell, temperature=300u"K")
ensemble = run_dynamics(sim, tspan, distribution; selection=1:nsamples,
    dt=0.1u"fs", output=OutputPosition, trajectories=nsamples, callback=terminate)
```

## Visualisation

To show the effect of the truncation procedure, we have run 20 trajectories with and without the truncation
function starting with an initial translation energy at 1.0 eV. For both figures, the total and kinetic energies are shown in
the top panels along with the H-H distance and centre of mass z coordinate for each
individual trajectory.

![truncation](../assets/figures/scattering_truncation.png)



```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
