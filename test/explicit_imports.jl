using ExplicitImports
using DiffEqPhysics
using Test

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(DiffEqPhysics) === nothing
    @test check_no_stale_explicit_imports(DiffEqPhysics) === nothing
end
