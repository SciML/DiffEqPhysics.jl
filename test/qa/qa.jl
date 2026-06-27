using SciMLTesting, DiffEqPhysics, Test

run_qa(
    DiffEqPhysics;
    explicit_imports = true,
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (
                :GradientConfig, :derivative, :gradient, :gradient!,  # ForwardDiff: not public
                :plot,                                                # RecipesBase: not public
            ),
        ),
    ),
)

# JET reports genuine errors in src/plot.jl (RecipesBase.plot has no inferable
# method for OrbitPlot, and DiffEqPhysics.plot / plot! are undefined in
# plot_orbits). The finding is JET-analysis-version dependent (JET 0.11 on the
# `1` lane surfaces it; JET 0.9, the only version installable on the `lts` lane,
# does not), so it is kept as a static tracked @test_broken rather than run_qa's
# version-auto-flagging jet_broken (which would Unexpected-Pass on lts and hard-FAIL
# on `1`). Pending fix, see https://github.com/SciML/DiffEqPhysics.jl/issues/118
@testset "JET (broken)" begin
    @test_broken false  # JET: no method `plot(::OrbitPlot)` / `DiffEqPhysics.plot` undefined in src/plot.jl — https://github.com/SciML/DiffEqPhysics.jl/issues/118
end
