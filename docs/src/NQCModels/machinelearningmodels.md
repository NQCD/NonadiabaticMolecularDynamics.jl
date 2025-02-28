```@setup logging
@info "Expanding src/NQCModels/machinelearningmodels.md..."
start_time = time()
```
# [Machine learning interatomic potentials](@id ml-pes-models)

## Loading ML models with `ase` calculators

Using the ASE interface within (NQCModels.jl)[https://github.com/NQCD/NQCModels.jl] we can directly use models trained e.g. using [MACE](https://github.com/ACEsuit/mace).

To use a MACE model, please load any pre-trained model into a given path you can access.

First we load the model into an `ase` calculator and attach it to our diatomic
hydrogen molecule.
```julia
using PythonCall

ase = pyimport("ase")
mace_calc = pyimport("mace.calculators")

h2 = ase.Atoms("H2", [(0, 0, 0), (0, 0, 0.74)])

calculator = mace_calc.MACECalculator(model_path="../assets/mace/h2cu.model", device="cpu", default_dtype="float32")
h2.set_calculator(calculator)
```

We can obtain the energies and forces from `ase` directly in the usual way, converting
them to atomic units using [UnitfulAtomic](https://github.com/sostock/UnitfulAtomic.jl).
```julia-repl
using Unitful, UnitfulAtomic;
austrip(h2.get_total_energy() * u"eV")
austrip.(h2.get_forces() .* u"eV/Å")
```

!!! warning

    Note that this is an arbitrary model not trained on H2, hence the calculation of the
    potential energy and forces most likely do not make sense.

Then, we can convert the ASE output into the format used in NQCModels,
which makes it possible to use the MACE model e.g. for molecular dynamics calculations
within NQCDynamics.jl:
```julia-repl
using NQCModels;
model = AdiabaticASEModel(h2);

r = [0 0; 0 0; 0 ustrip(auconvert(0.74u"Å"))]

potential(model, r)
derivative(model, r)
```

## Example: Loading SchNet (SchNetPack) models through `ase`

Using the [ASE interface](@ref) we can directly use models trained using
[SchNetPack](https://github.com/atomistic-machine-learning/schnetpack).

!!! warning

    The examples on this page do not run during the documentation build due to `schnetpack`
    causing segfaults when installed in the build environment.
    The causes of this is not currently clear but we have temporarily disabled these examples
    in the build.

    However, the examples should still be correct and you are welcome to try them with
    your own schnetpack trained models.

To use a SchNet model, please load any pre-trained model into a given path you can access.
Here, our SchNet model is named "best_model" as is common in SchNet and
provide the relative path.

First we load the model into an `ase` calculator and attach it to our diatomic
hydrogen molecule.

```julia
using PythonCall

ase = pyimport("ase")
spkutils = pyimport("schnetpack.utils")
spkinterfaces = pyimport("schnetpack.interfaces")

spk_model = spkutils.load_model("../assets/schnetpack/best_model"; map_location="cpu")

h2 = ase.Atoms("H2", [(0, 0, 0), (0, 0, 0.74)])

calc = spkinterfaces.SpkCalculator(spk_model, energy="energy", forces="forces")
h2.set_calculator(calc)
```

We can obtain the energies and forces from `ase` directly in the usual way, converting
them to atomic units using [UnitfulAtomic](https://github.com/sostock/UnitfulAtomic.jl).


```julia-repl
using Unitful, UnitfulAtomic;
austrip(pyconvert(Float64,h2.get_total_energy()) * u"eV")
austrip.(pyconvert(Matrix{Float64}, h2.get_forces()) .* u"eV/Å")
austrip(pyconvert(Float64,h2.get_total_energy()) * u"eV")
austrip.(pyconvert(Matrix{Float64}, h2.get_forces()) .* u"eV/Å")
```

!!! warning

    Note that this is an arbitrary model not trained on H2, hence the calculation of the
    potential energy and forces most likely do not make sense.

Then, we can convert the ASE output into the format used in NQCModels,
which makes it possible to use the SchNet model e.g. for molecular dynamics calculations
within NQCDynamics.jl:

```julia-repl
using NQCModels;
model = AdiabaticASEModel(h2);

r = [0 0; 0 0; 0 ustrip(auconvert(0.74u"Å"))]

potential(model, r)
derivative(model, r)
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
