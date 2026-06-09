using DiffEqPhysics, Aqua, JET, Test

@testset "Aqua" begin
    Aqua.test_all(DiffEqPhysics)
end

@testset "JET" begin
    JET.test_package(DiffEqPhysics; target_defined_modules = true)
end
