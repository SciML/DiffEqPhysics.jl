using SafeTestsets, Test, Pkg

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Core"
    @safetestset "Hamiltonian Test" begin
        include("hamiltonian_test.jl")
    end
    #@safetestset "N-Body Test" begin include("nbody_test.jl") end

    @safetestset "Explicit Imports" begin
        include("explicit_imports.jl")
    end
end

if GROUP == "QA"
    Pkg.activate("qa")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    include("qa/qa.jl")
end
