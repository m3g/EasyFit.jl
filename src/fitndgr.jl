#
# Ndgr fit: fit to an n-th degree polynomial
#

struct Ndgr{T} <: Fit{T}
    lscoeff::Vector{T}
    R::T
    x::Vector{T}
    y::Vector{T}
    ypred::Vector{T}
    residues::Vector{T}
end

"""
    fitndgr(x,y,n)

Obtains the quadratic fit: ``y = p[n+1]*x^n + p[n]*x^(n-1) + ... + p[2]*x + p[1]``

Optional lower and upper bounds for p[i] can be provided using two arrays of length n+1, for example:

```
fitndgr(x,y,4; l=[-1,-1.,-1,-1], u=[5.,7.,8.,7.))
```
"""
function fitndgr(
    X::AbstractArray{T1}, Y::AbstractArray{T2}, n::T3;
    l=[-Inf], u=[Inf], c=nothing,
    options::Options=Options()
) where {T1<:Real, T2<:Real, T3<:Int}
    # Check data
    X, Y, data_type = checkdata(X, Y, options)
    # Set bounds
    if length(l) == n+1
        lower = l
    else
        lower = fill(l[1],n+1)
    end
    lower = convert(typeof(Y), lower)
    if length(u) == n+1
        upper = u
    else
        upper = fill(u[1],n+1)
    end
    upper = convert(typeof(Y), upper)
    if isnothing(c)
        # Set model
        model(x, p) = sum([p[i] .* (x.^(n+1-i)) for i in 1:n+1])
        # Fit
        fit = find_best_fit(model, X, Y, n+1, options, lower, upper)
        # Analyze results and return
        R = pearson(X, Y, model, fit)
        x, y, ypred = finexy(X, options.fine, model, fit)
        return Ndgr(fit.param, R, x, y, ypred, fit.resid)
    else
        lower = lower[1:n]
        upper = upper[1:n]
        # Set model
        model_const(x, p) = sum([p[i] .* (x.^(n+1-i)) for i in 1:n]) + c
        # Fit
        fit = find_best_fit(model_const, X, Y, n, options, lower, upper)
        # Analyze results and return
        R = pearson(X, Y, model_const, fit)
        x, y, ypred = finexy(X, options.fine, model_const, fit)
        return Ndgr([fit.param..., c], R, x, y, ypred, fit.resid)
    end
end

function (fit::EasyFit.Ndgr)(x::Real)
    n = length(fit.lscoeff)-1
    return sum([fit.lscoeff[i] * x^(n+1-i) for i in 1:n+1])
end

function Base.show(io::IO, fit::Ndgr)
    println(io,
        """
        ------------------- n-th degree Fit -----------------

        Equation: y = sum([p[i] * x^(n+1-i) for i in 1:n+1])

        With: p = $(fit.lscoeff)

        Pearson correlation coefficient, R = $(fit.R))
        Average square residue = $(mean(fit.residues .^ 2))

        Predicted Y: ypred = [ $(fit.ypred[1]), $(fit.ypred[2]), ...]
        residues = [ $(fit.residues[1]), $(fit.residues[2]), ...]

        -----------------------------------------------"""
    )
end

export fitndgr

@testitem "fitndgr" begin
    using Statistics: mean
    x = sort(rand(10))
    y = @. 6 * x^4 + 4 * x^3 + 3 * x^2 + 2 * x + 1
    f = fitndgr(x, y, 4)
    @test f.R ≈ 1
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R^2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred == f.(x))

    x = Float32.(x)
    y = Float32.(y)
    f = fitndgr(x, y, 4)
    @test typeof(f.R) == Float32
end
