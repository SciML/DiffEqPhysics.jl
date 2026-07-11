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
    api_docs_kwargs = (;
        rendered = true,
        ignore = (
            :DEFAULT_PLOT_FUNC, :add_labels!, :plot_indices, :plottable_indices,
            :recursive_bottom_eltype, :recursive_mean, :solplot_vecs_and_labels,
            :tuples, :vecarr_to_vectors, :vecvec_to_mat,
        ),
        rendered_ignore = (
            :AP, :AbstractAnalyticalProblem, :AbstractDiffEqArray,
            :AbstractVectorOfArray, :AddVector, :AffineOperator, :AllObserved,
            :AnalyticalProblem, :ArrayPartition, :BVPFunction, :BVProblem,
            :BatchIntegralFunction, :BlockDiagonalOperator, :BrownFullBasicInit,
            :CallbackSet, :CheckInit, :Clocks, :ContinuousCallback, :DAEFunction,
            :DAEProblem, :DAESolution, :DDEFunction, :DDEProblem,
            :DEFAULT_PLOT_FUNC, :DEVerbosity, :DefaultInit, :DiagonalOperator,
            :DiffEqArray, :DiffEqBase, :DiscreteCallback, :DiscreteFunction,
            :DiscreteProblem, :DynamicalBVPFunction, :DynamicalDDEFunction,
            :DynamicalDDEProblem, :DynamicalODEFunction, :DynamicalODEProblem,
            :DynamicalSDEFunction, :DynamicalSDEProblem, :EigenvalueProblem,
            :EigenvalueSolution, :EigenvalueTarget, :EnsembleAnalysis,
            :EnsembleContext, :EnsembleDistributed, :EnsembleProblem,
            :EnsembleSerial, :EnsembleSolution, :EnsembleSplitThreads,
            :EnsembleSummary, :EnsembleTestSolution, :EnsembleThreads,
            :FunctionOperator, :HomotopyNonlinearFunction, :HomotopyProblem,
            :IdentityOperator, :ImplicitDiscreteFunction, :ImplicitDiscreteProblem,
            :IncrementingODEFunction, :IncrementingODEProblem, :IntegralFunction,
            :IntegralProblem, :IntegralSolution, :IntervalNonlinearFunction,
            :IntervalNonlinearProblem, :InvertibleOperator, :LinearAliasSpecifier,
            :LinearProblem, :LinearSolution, :MatrixOperator,
            :MultiObjectiveOptimizationFunction, :NamedArrayPartition, :NoInit,
            :NoiseProblem, :NonlinearFunction, :NonlinearLeastSquaresProblem,
            :NonlinearProblem, :NonlinearSolution, :NullOperator,
            :ODEAliasSpecifier, :ODEFunction, :ODEInputFunction, :ODEProblem,
            :ODESolution, :OptimizationFunction, :OptimizationProblem,
            :OptimizationSolution, :OverrideInit, :PDENoTimeSolution, :PDEProblem,
            :PDETimeSeriesSolution, :RODEFunction, :RODEProblem, :RODESolution,
            :RecursiveArrayTools, :ReturnCode, :SCCNonlinearProblem, :SDDEFunction,
            :SDDEProblem, :SDEFunction, :SDEProblem, :SampledIntegralProblem,
            :ScalarOperator, :SciMLBase, :SciMLOperators, :SecondOrderBVProblem,
            :SecondOrderDDEProblem, :SecondOrderODEProblem,
            :SensitivityADPassThrough, :ShampineCollocationInit, :SplitFunction,
            :SplitODEProblem, :SplitSDEFunction, :SplitSDEProblem, :StaticWOperator,
            :SteadyStateProblem, :SteadyStateSolution, :TensorProductOperator,
            :TensorSumOperator, :TimeDomain, :TwoPointBVPFunction, :TwoPointBVProblem,
            :TwoPointDynamicalBVPFunction, :TwoPointSecondOrderBVProblem, :VA,
            :VectorContinuousCallback, :VectorOfArray, :WOperator, :add_labels!,
            :add_saveat!, :add_tstop!, :addat!, :addat_non_user_cache!, :addsteps!,
            :auto_dt_reset!, :cache_operator, :change_t_via_interpolation!,
            :check_error, :check_keywords, :concretize, :copyat_or_push!,
            :deleteat!, :deleteat_non_user_cache!, :derivative_discontinuity!,
            :diffeq_to_arrays, :discretize, :du_cache, :finalize!, :first_tstop,
            :full_cache, :get_dt, :get_du, :get_du!, :get_proposed_dt, :get_rng,
            :get_tmp_cache, :getindepsym_defaultt, :has_adjoint, :has_concretization,
            :has_exp, :has_expmv, :has_expmv!, :has_ldiv, :has_ldiv!, :has_mul,
            :has_mul!, :has_rng, :has_tstop, :init, :initialize!, :interpret_vars,
            :is_discrete_time_domain, :iscached, :isclock, :isconstant, :iscontinuous,
            :isconvertible, :isdiscrete, :isinplace, :islinear, :issolverstepclock,
            :issquare, :kronsum, :plot_indices, :plottable_indices, :pop_tstop!,
            :rand_cache, :ratenoise_cache, :recursive_bottom_eltype,
            :recursive_mean, :recursive_one, :recursive_unitless_bottom_eltype,
            :recursive_unitless_eltype, :recursivecopy, :recursivecopy!,
            :recursivecopyto!, :recursivefill!, :reeval_internals_due_to_modification!,
            :reinit!, :remake, :resize_non_user_cache!, :savevalues!, :set_abstol!,
            :set_proposed_dt!, :set_reltol!, :set_rng!, :set_t!, :set_u!,
            :solplot_vecs_and_labels, :solve, :solve!, :step!, :supports_solve_rng,
            :symbolic_discretize, :terminate!, :tuples, :u_cache, :u_modified!,
            :update_coefficients, :update_coefficients!, :user_cache, :vecarr_to_vectors,
            :vecvec_to_mat, :vecvecapply, :warn_compat,
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
