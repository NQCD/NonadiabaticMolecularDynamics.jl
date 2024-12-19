```@setup logging
@info "Expanding src/NQCModels/ase.md..."
start_time = time()
```
# ASE interface

The easiest way to obtain potentials and forces from established codes is to
use the interfaces implemented in [ASE](https://wiki.fysik.dtu.dk/ase/).

We provide the [`AdiabaticASEModel`](@ref) which wraps an ASE atoms object and its
associated calculator to implement the required [`potential`](@ref) and
[`derivative`](@ref) functions.

!!! note

    The interface works by calling the relevant Python functions using
    [PyCall](https://github.com/JuliaPy/PyCall.jl).
    To use PyCall, you must make sure that your `python` version contains all the
    relevant packages, such as ase.
    [PyCall](https://github.com/JuliaPy/PyCall.jl) can be configured to use a particular
    pre-installed Python or install its own.
    Refer to the [PyCall README](https://github.com/JuliaPy/PyCall.jl) for installation
    and configuration instructions.
    
## Example

First, it is necessary to import `ase` and create the `ase.Atoms` object and attach
the desired calculator. This works exactly as in Python:
```@example ase
using PythonCall: pyimport, pylist
using ASEconvert

emt = pyimport("ase.calculators.emt")

h2 = ase.Atoms("H2", pylist([(0, 0, 0), (0, 0, 0.74)]))
h2.calc = emt.EMT()
nothing # hide
```

Next, the [`AdiabaticASEModel`](@ref) is created by passing the `ase.Atoms` object directly
to the model:
```@repl ase
using NQCModels
model = AdiabaticASEModel(h2)
```
Now the model can be used in the same way as any of the previously introduced
analytic models.
```@repl ase
potential(model, rand(3, 2))
derivative(model, rand(3, 2))
```

!!! tip 

    In theory, this should work with any of the ASE calculators that correctly implement
    the `get_potential_energy` and `get_forces` functions. For instance, you can use
    [SchNetPack (SPK)](https://github.com/atomistic-machine-learning/schnetpack) by
    passing their ASE calculator to the `AdiabaticASEModel`.
    Take a look at [machine learning models](@ref ml-pes-models) to learn more.
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
