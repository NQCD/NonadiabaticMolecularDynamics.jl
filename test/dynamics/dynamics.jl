using SafeTestsets
using Test

@time @safetestset "Classical Tests" begin include("classical.jl") end
@time @safetestset "Langevin Tests" begin include("langevin.jl") end
@time @safetestset "MDEF BAOAB Tests" begin include("mdef_baoab.jl") end
@time @safetestset "MDEF Tests" begin include("mdef.jl") end
@time @safetestset "RPMDEF Tests" begin include("rpmdef.jl") end
@time @safetestset "BCBwithTsit5 Tests" begin include("bcbwithtsit5.jl") end
@time @safetestset "FSSH Tests" begin include("fssh.jl") end
@time @safetestset "NRPMD Tests" begin include("nrpmd.jl") end
@time @safetestset "CMM Tests" begin include("cmm.jl") end
@time @safetestset "Cell Boundary Callback Tests" begin include("cell_boundary_callback.jl") end
@time @safetestset "Ehrenfest Tests" begin include("ehrenfest.jl") end
@time @safetestset "IESH Tests" begin include("iesh.jl") end
