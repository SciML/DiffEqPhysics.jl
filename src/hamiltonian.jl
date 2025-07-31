const TOO_MANY_ARGUMENTS_ERROR_MESSAGE = """
                                         All methods for the model function `H` had too many arguments. A
                                         Hamiltonian `H` must define either `H(p, q, param)` or `H(p, q, param, t)`. This error
                                         can be thrown if you define an Hamiltonian for example as `H(p, q, param1, param2, t)`.
                                         For more information on the required number of arguments for the function
                                         you were defining, consult the documentation for the `HamiltonianProblem`.
                                         """

struct HamiltonainTooManyArgumentsError <: Exception
    fname::String
    f::Any
end

function Base.showerror(io::IO, e::HamiltonainTooManyArgumentsError)
    println(io, TOO_MANY_ARGUMENTS_ERROR_MESSAGE)
    print(io, "Offending function: ")
    printstyled(io, e.fname; bold = true, color = :red)
    println(io, "\nMethods:")
    println(io, methods(e.f))
end

const TOO_FEW_ARGUMENTS_ERROR_MESSAGE = """
                                        All methods for the Hamiltonian `H` had too few arguments. A
                                        Hamiltonian `H` must define either `H(p, q, param)` or `H(p, q, param, t)`. This error
                                        can be thrown if you define an Hamiltonian for example as `H(p, q)`.
                                        Note that `param` must be in the arguments list even if it's not used.
                                        For more information on the required number of arguments for the function
                                        you were defining, consult the documentation for the `HamiltonianProblem`.
                                        """

struct HamiltonainTooFewArgumentsError <: Exception
    fname::String
    f::Any
end

function Base.showerror(io::IO, e::HamiltonainTooFewArgumentsError)
    println(io, TOO_FEW_ARGUMENTS_ERROR_MESSAGE)
    print(io, "Offending function: ")
    printstyled(io, e.fname; bold = true, color = :red)
    println(io, "\nMethods:")
    println(io, methods(e.f))
end

const ARGUMENTS_ERROR_MESSAGE = """
                                Methods dispatches for the Hamiltonian `H` do not match the required number.
                                A Hamiltonian `H` must define either `H(p, q, param)` or `H(p, q, param, t)`. This error
                                can be thrown if you define an Hamiltonian for example as `H(p)`.
                                Note that arguments must be in the arguments list even if it's not used.
                                For more information on the required number of arguments for the function
                                you were defining, consult the documentation for the `HamiltonianProblem`.
                                """

struct HamiltonainFunctionArgumentsError <: Exception
    fname::String
    f::Any
end

function Base.showerror(io::IO, e::HamiltonainFunctionArgumentsError)
    println(io, ARGUMENTS_ERROR_MESSAGE)
    print(io, "Offending function: ")
    printstyled(io, e.fname; bold = true, color = :red)
    println(io, "\nMethods:")
    println(io, methods(e.f))
end

struct HamiltonianProblem{iip} <: DiffEqBase.AbstractDynamicalODEProblem end

"""
    HamiltonianProblem(H, p0, q0, tspan, param=nothing; kwargs...)
    HamiltonianProblem((dp, dq), p0, q0, tspan, param=nothing; kwargs...)

Define a physical system by its Hamiltonian function `H(p, q, param)` or the function pair
`dp = -∂H/∂q` and `dq = ∂H/∂p`.

The equations of motion are then given by `q̇ = ∂H/∂p = dq` and `ṗ  = -∂H/∂q = dp`.

The initial values for canonical impulses `p0` and coordinates `q0` may be scalars, `SArray`s or
other `AbstractArray`s. Their type determines the type of functions `dp` and `dq`.

In the first two cases the derivative functions require the signatures
`dp(p, q, param, t)` and `dq(p, q, param, t)` while
in the latter case the partial derivatives use mutating signatures
`dp!(Δp, p, q, param, t)` and `dq!(Δq, p, q, param, t)` with predefined arrays `Δp` and `Δq`.

If the Hamiltonian function is given, `dp` and `dq` are calculated automatically using
AD (`ForwardDiff`).

!!! note


`H` may be defined with or without time as fourth argument. If both methods are defined,
that with 4 arguments is used.
"""
function HamiltonianProblem(
        H, p0::S, q0::T, tspan, param = NullParameters(); kwargs...) where {S, T}
    iip = T <: AbstractArray && !(T <: StaticArraysCore.SArray) && S <: AbstractArray &&
          !(S <: StaticArraysCore.SArray)
    HamiltonianProblem{iip}(H, p0, q0, tspan, param; kwargs...)
end

struct PhysicsTag end

function generic_derivative(q0, hami, x)
    ForwardDiff.gradient(hami, x)
end

function generic_derivative(q0::Number, hami, x)
    ForwardDiff.derivative(hami, x)
end

function HamiltonianProblem{false}(
        (dp, dq)::Tuple{Any, Any}, p0, q0, tspan, param = NullParameters(); kwargs...)
    return ODEProblem(DynamicalODEFunction{false}(dp, dq),
        ArrayPartition(p0, q0), tspan, param; kwargs...)
end
function HamiltonianProblem{true}(
        (dp, dq)::Tuple{Any, Any}, p0, q0, tspan, param = NullParameters(); kwargs...)
    return ODEProblem(
        DynamicalODEFunction{true}(dp, dq), ArrayPartition(p0, q0), tspan, param; kwargs...)
end

function HamiltonianProblem{false}(H, p0, q0, tspan, param = NullParameters(); kwargs...)
    try
        isinplace(H, 4)
    catch e
        if e isa SciMLBase.TooManyArgumentsError
            throw(HamiltonainTooManyArgumentsError(e.fname, e.f))
        elseif e isa SciMLBase.TooFewArgumentsError
            throw(HamiltonainTooFewArgumentsError(e.fname, e.f))
        elseif e isa SciMLBase.FunctionArgumentsError
            throw(HamiltonainFunctionArgumentsError(e.fname, e.f))
        end
    end

    if 4 in DiffEqBase.numargs(H)
        dp = (p, q, param, t) -> generic_derivative(q0, q -> -H(p, q, param, t), q)
        dq = (p, q, param, t) -> generic_derivative(q0, p -> H(p, q, param, t), p)
    else
        issue_depwarn()
        dp = (p, q, param, t) -> generic_derivative(q0, q -> -H(p, q, param), q)
        dq = (p, q, param, t) -> generic_derivative(q0, p -> H(p, q, param), p)
    end
    return HamiltonianProblem{false}((dp, dq), p0, q0, tspan, param; kwargs...)
end

function HamiltonianProblem{true}(H, p0, q0, tspan, param = NullParameters(); kwargs...)
    try
        isinplace(H, 4)
    catch e
        if e isa SciMLBase.TooManyArgumentsError
            throw(HamiltonainTooManyArgumentsError(e.fname, e.f))
        elseif e isa SciMLBase.TooFewArgumentsError
            throw(HamiltonainTooFewArgumentsError(e.fname, e.f))
        elseif e isa SciMLBase.FunctionArgumentsError
            throw(HamiltonainFunctionArgumentsError(e.fname, e.f))
        end
    end
    let cp = ForwardDiff.GradientConfig(PhysicsTag(), p0),
        cq = ForwardDiff.GradientConfig(PhysicsTag(), q0), vfalse = Val(false)

        if 4 in DiffEqBase.numargs(H)
            dp = (Δp, p, q, param,
                t) -> ForwardDiff.gradient!(Δp, q->-H(p, q, param, t), q, cq, vfalse)
            dq = (Δq, p, q, param,
                t) -> ForwardDiff.gradient!(Δq, p -> H(p, q, param, t), p, cp, vfalse)
        else
            issue_depwarn()
            dp = (Δp, p, q, param,
                t) -> ForwardDiff.gradient!(Δp, q->-H(p, q, param), q, cq, vfalse)
            dq = (Δq, p, q, param,
                t) -> ForwardDiff.gradient!(Δq, p -> H(p, q, param), p, cp, vfalse)
        end
        return HamiltonianProblem{true}((dp, dq), p0, q0, tspan, param; kwargs...)
    end
end

function issue_depwarn()
    Base.depwarn(
        "Hamiltonians with 3 arguments are deprecated; please use `H(p, q, params, t)`",
        :HamiltonianProblem)
end
