```@setup logging
@info "Expanding src/api/NQCModels/FrictionProviders.md..."
start_time = time()
```

# FrictionProviders.jl

This package provides MLIP interfaces for the prediction of electronic friction tensors for molecular systems. 

```@autodocs
Modules=[FrictionProviders]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
