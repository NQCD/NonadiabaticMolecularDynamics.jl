```@setup logging
@info "Expanding src/api/NQCModels/mace.md..."
start_time = time()
```

# MACEModels.jl

```@autodocs
Modules = [MACEModels, MACEModels.Ensemble]
Private = true
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
