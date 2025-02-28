```@setup logging
@info "Expanding src/api/NQCModels/frictionmodels.md..."
start_time = time()
```

# FrictionModels

```@autodocs
Modules=[NQCModels.FrictionModels]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
