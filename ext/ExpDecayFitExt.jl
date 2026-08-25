# requires JuMP and Ipopt
module ExpDecayFitExt

using TestItems
using JuMP
import Ipopt
using Statistics: mean, cor
using Unitful: ustrip
import EasyFit: Fit, Options, fitexpdecay

#
# Normalized multiple-exponential decay
#
# y(t) = sum(a[i] * exp(-t/b[i]) for i in 1:n) + c
#
# subject to sum(a) + c == 1 (so that y(0) == 1) and b[i] > 0 for all i. The
# constant `c` can optionally be fixed to a user-provided value, otherwise it
# is fitted freely subject to c >= 0.
#
struct ExpDecayFit{T} <: Fit{T}
    n::Int
    a::Vector{T}
    b::Vector{T}
    c::T
    R2::T
    x::Vector{T}
    y::Vector{T}
    ypred::Vector{T}
    residues::Vector{T}
end

normexp_model(t::AbstractVector{<:Real}, a::AbstractVector{<:Real}, b::AbstractVector{<:Real}, c::Real) =
    [sum(a[j] * exp(-ti / b[j]) for j in eachindex(a)) + c for ti in t]

# Converts any AbstractVector (including OffsetArrays) into a standard 1-based Vector{Float64}.
# Note: a comprehension over `eachindex(v)` would inherit v's (possibly offset) axes,
# so the output vector is built by explicit position instead.
function to_data_vector(v::AbstractVector)
    out = Vector{Float64}(undef, length(v))
    for (k, i) in enumerate(eachindex(v))
        out[k] = Float64(ustrip(v[i]))
    end
    return out
end

# By default the time of each data point is its index step, i.e. the first
# data point is assumed to correspond to t = 0.
default_times(Y::AbstractVector) = collect(Float64, 0:(length(Y)-1))

function solve_decay_nlp(t::Vector{Float64}, y::Vector{Float64}, n::Int, c_fixed::Union{Nothing,Real}, options::Options)
    tmin, tmax = extrema(t)
    trange = max(tmax - tmin, 1.0)
    bmin = 1e-6 * trange

    # The model is built once and re-solved from different starting points
    # (changing only the variables' start values) instead of being rebuilt
    # on every multistart trial.
    jumpmodel = Model(Ipopt.Optimizer)
    set_silent(jumpmodel)
    @variable(jumpmodel, a[1:n])
    @variable(jumpmodel, b[1:n] >= bmin)
    if isnothing(c_fixed)
        @variable(jumpmodel, c >= 0)
        @constraint(jumpmodel, sum(a) + c == 1)
        @NLobjective(jumpmodel, Min,
            sum((y[k] - (sum(a[j] * exp(-t[k] / b[j]) for j in 1:n) + c))^2 for k in eachindex(t)))
    else
        @constraint(jumpmodel, sum(a) == 1 - c_fixed)
        @NLobjective(jumpmodel, Min,
            sum((y[k] - (sum(a[j] * exp(-t[k] / b[j]) for j in 1:n) + c_fixed))^2 for k in eachindex(t)))
    end

    best_obj = +Inf
    local best_a::Vector{Float64}, best_b::Vector{Float64}, best_c::Float64
    nbest = 0
    ntrial = 0
    while nbest < options.nbest && ntrial < options.maxtrials
        ntrial += 1
        c0 = isnothing(c_fixed) ? 0.0 : c_fixed
        a0 = rand(n)
        a0 .*= (1 - c0) / sum(a0)
        b0 = trange .* (10.0 .^ (4 .* rand(n) .- 2))
        set_start_value.(a, a0)
        set_start_value.(b, b0)
        isnothing(c_fixed) && set_start_value(c, c0)
        try
            optimize!(jumpmodel)
            status_ok = termination_status(jumpmodel) in (JuMP.MOI.LOCALLY_SOLVED, JuMP.MOI.OPTIMAL) &&
                        primal_status(jumpmodel) == JuMP.MOI.FEASIBLE_POINT
            if status_ok
                obj = objective_value(jumpmodel)
                improve = obj - best_obj
                if improve < options.besttol
                    if improve < 0
                        nbest = 1
                        best_obj = obj
                        best_a = value.(a)
                        best_b = value.(b)
                        best_c = isnothing(c_fixed) ? value(c) : Float64(c_fixed)
                    else
                        nbest += 1
                    end
                end
            end
        catch err
            options.debug && @warn "Ipopt trial failed" exception = err
        end
    end
    if nbest == 0
        error("""
        Could not obtain any successful fit, probably the data is not well posed.
        Further information can be obtained by adding `options=Options(debug=true)` as kwarg.
        """)
    end
    return best_a, best_b, best_c
end

"""
    fitexpdecay(Y; n::Int=1, t=nothing, c=nothing, options::Options=Options())

Fits a normalized multiple-exponential decay model to `Y`:

``y(t) = \\sum_{i=1}^n a_i \\, e^{-t/b_i} + c``

subject to the constraints ``\\sum_i a_i + c = 1`` (so that ``y(0) = 1``) and
``b_i > 0`` for all decay rates. The fit is solved as a constrained nonlinear
least-squares problem using JuMP with the Ipopt solver, requiring
`using JuMP, Ipopt` to be loaded.

`Y` can be a plain vector or an `OffsetArray` (for example with axis
`0:length(Y)-1`). By default the time associated to each data point is
`index - first(index)`, i.e. the first data point is assumed to correspond
to `t = 0`. A time vector `t` of the same length as `Y` can optionally be
provided explicitly.

The independent constant `c` is fitted freely (subject to `c >= 0`, and hence
`sum(a) <= 1`) by default. It can optionally be fixed to a user-provided value
with the `c` keyword, e.g. `c=0.0`; a fixed value is not required to be
non-negative (`sum(a)` is then set to `1 - c` accordingly).

# Examples

```jldoctest
julia> using JuMP, Ipopt

julia> t = 0:0.1:5; y = @. 0.7*exp(-t/0.5) + 0.3*exp(-t/3);

julia> fit = fitexpdecay(y; t=collect(t), n=2)
```
"""
function fitexpdecay(
    Y::AbstractVector{<:Real};
    n::Int=1,
    t::Union{Nothing,AbstractVector{<:Real}}=nothing,
    c::Union{Nothing,Real}=nothing,
    options::Options=Options(),
)
    y = to_data_vector(Y)
    tv = isnothing(t) ? default_times(Y) : to_data_vector(t)
    length(tv) == length(y) || throw(ArgumentError("t and Y must have the same length."))
    a, b, cfit = solve_decay_nlp(tv, y, n, c, options)
    ind = sortperm(b)
    a = a[ind]
    b = b[ind]
    ypred = normexp_model(tv, a, b, cfit)
    residues = ypred .- y
    R = cor(y, ypred)^2
    tmin, tmax = extrema(tv)
    xfine = collect(range(tmin, tmax, length=options.fine))
    yfine = normexp_model(xfine, a, b, cfit)
    return ExpDecayFit(n, a, b, cfit, R, xfine, yfine, ypred, residues)
end

function (fit::ExpDecayFit)(x::Real)
    return sum(fit.a[i] * exp(-x / fit.b[i]) for i in eachindex(fit.a)) + fit.c
end

function Base.show(io::IO, fit::ExpDecayFit)
    println(io,
        """
        -------- Normalized multiple-exponential decay fit --------

        Equation: y = sum(a[i] exp(-t/b[i]) for i in 1:$(fit.n)) + c, with sum(a) + c = 1 and b .> 0

        With: a = $(fit.a)
              b = $(fit.b)
              c = $(fit.c)

        Correlation coefficient, R² = $(fit.R2)
        Average square residue = $(mean(fit.residues .^ 2))

        Predicted Y: ypred = [$(fit.ypred[1]), $(fit.ypred[2]), ...]
        residues = [$(fit.residues[1]), $(fit.residues[2]), ...]

        -------------------------------------------------------------"""
    )
end

@testitem "fitexpdecay" begin
    using JuMP, Ipopt
    using OffsetArrays
    using Statistics: mean

    t = 0:0.1:5
    y = @. 0.7 * exp(-t / 0.5) + 0.3 * exp(-t / 3)

    fit = fitexpdecay(y; t=collect(t), n=2)
    @test fit.R2 > 0.99
    @test isapprox(sum(fit.a) + fit.c, 1.0, atol=1e-6)
    @test all(fit.b .> 0)
    @test fit.c >= -1e-6
    @test all(fit.ypred - y .== fit.residues)
    @test all(isapprox.(fit.ypred, fit.(t), atol=1e-6))
    @test isapprox(fit(0.0), 1.0, atol=1e-6) # y(0) == 1 by construction

    # fixed constant: y(0) == 1 always holds by construction, so data with a
    # genuine baseline must itself be normalized that way (sum(a) is then 1 - c)
    y2 = @. 0.4 + 0.6 * exp(-t / 1.5)
    fit2 = fitexpdecay(y2; t=collect(t), n=1, c=0.4)
    @test fit2.c == 0.4
    @test isapprox(sum(fit2.a), 0.6, atol=1e-3)
    @test fit2.R2 > 0.99
    @test isapprox(fit2(0.0), 1.0, atol=1e-6)

    # fixing c to an arbitrary (even negative) value must not error; sum(a) is
    # then forced to 1 - c so that y(0) == 1 still holds regardless of c
    fit2alt = fitexpdecay(y; t=collect(t), n=2, c=-1.0)
    @test fit2alt.c == -1.0
    @test isapprox(sum(fit2alt.a), 2.0, atol=1e-6)

    # unconstrained least-squares would drive c negative here; the free fit must clamp c >= 0
    y3 = y .- 0.05
    fit3 = fitexpdecay(y3; t=collect(t), n=2)
    @test fit3.c >= -1e-6

    # default (index-based) times, including OffsetArray input
    # note: fitdefault and fitoffset are independent multistart optimizations
    # (random initial guesses), so they only agree up to solver tolerance,
    # not machine precision.
    yoff = OffsetArray(collect(y), 0:length(y)-1)
    fitdefault = fitexpdecay(collect(y); n=2)
    fitoffset = fitexpdecay(yoff; n=2)
    @test isapprox(fitdefault.ypred, fitoffset.ypred, atol=1e-3)
end

end # module
