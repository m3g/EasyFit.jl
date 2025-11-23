#
# Ndgr fit: fit to an n-th degree polynomial
#
struct Ndgr{T} <: Fit{T}
    lscoeff::Vector{T}
    R2::T
    x::Vector{T}
    y::Vector{T}
    ypred::Vector{T}
    residues::Vector{T}
end

"""
    fitndgr(x,y,n)

Obtains the polynomial fit of degree `n`: `y = p[n+1]*x^n + p[n]*x^(n-1) + ... + p[2]*x + p[1]`

Optional lower and upper bounds for p[i] can be provided using two arrays of length n+1, for example:

```
fitndgr(x,y,4; l=fill(-1.0, 5), u=[5.0, 7.0, 8.0, 7.0, 5.0])
```

"""
function fitndgr(
    X::AbstractArray{T1}, Y::AbstractArray{T2}, n::T3;
    l=nothing, u=nothing, c=nothing,
    options::Options=Options()
) where {T1<:Real,T2<:Real,T3<:Int}
    if n > 3
        # Check data
        X, Y, data_type = checkdata(X, Y, options)
        # Set bounds
        if isnothing(l)
            lower = fill(-Inf, n + 1)
        elseif length(l) == n + 1
            lower = l
        else
            throw(ArgumentError("Length of lower bound l must be n + 1, got $(length(l))"))
        end
        if isnothing(u)
            upper = fill(+Inf, n + 1)
        elseif length(u) == n + 1
            upper = u
        else
            throw(ArgumentError("Length of upper bound u must be n + 1, got $(length(u))"))
        end
        lower = convert(typeof(Y), lower)
        upper = convert(typeof(Y), upper)
        if isnothing(c)
            # Set model
            model(x, p) = sum(p[i] * (x .^ (i - 1)) for i in n+1:-1:1)
            # Fit
            fit = find_best_fit(model, X, Y, n + 1, options, lower, upper)
            # Analyze results and return
            R = R2(X, Y, model, fit)
            x, y, ypred = finexy(X, options.fine, model, fit)
            return Ndgr(fit.param, R, x, y, ypred, fit.resid)
        else
            lower = lower[1:n]
            upper = upper[1:n]
            # Set model
            model_const(x, p) = sum(p[i] * (x^i) for i in n:-1:1) + c
            # Fit
            fit = find_best_fit(model_const, X, Y, n, options, lower, upper)
            # Analyze results and return
            R = R2(X, Y, model_const, fit)
            x, y, ypred = finexy(X, options.fine, model_const, fit)
            return Ndgr([fit.param..., c], R, x, y, ypred, fit.resid)
        end
    else
        if isnothing(l)
            l = EasyFit.lower()
        end
        if isnothing(u)
            u = EasyFit.upper()
        end
        if n == 1
            return fitlinear(
                X, Y;
                l=l, u=u, b=c,
                options=options
            )
        elseif n == 2
            return fitquadratic(
                X, Y;
                l=l, u=u, c=c,
                options=options
            )
        elseif n == 3
            return fitcubic(
                X, Y;
                l=l, u=u, d=c,
                options=options
            )
        end
    end
end

function (fit::EasyFit.Ndgr)(x::Real)
    n = length(fit.lscoeff) - 1
    return sum(fit.lscoeff[i] * x^(i - 1) for i in n+1:-1:1)
end

function Base.show(io::IO, fit::Ndgr)
    println(io, chomp(
        """
        ------------- n-th degree polynomial degree Fit -------------

        Equation: y = sum(p[i] * x^(i-1) for i in n+1:-1:1)

        With: p = $(fit.lscoeff)

        Correlation coefficient, R² = $(fit.R2)
        Average square residue = $(mean(fit.residues .^ 2))

        Predicted Y: ypred = [ $(fit.ypred[1]), $(fit.ypred[2]), ...]
        residues = [ $(fit.residues[1]), $(fit.residues[2]), ...]

        -------------------------------------------------------------
        """)
    )
end

export fitndgr

@testitem "fitndgr for n = 1" begin
    using Statistics: mean
    x = sort(rand(10))
    y = @. 2 * x + 1
    f = fitndgr(x, y, 1)
    @test f.R2 ≈ 1
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred == f.(x))

    x = Float32.(x)
    y = Float32.(y)
    f = fitndgr(x, y, 1)
    @test typeof(f.R2) == Float32
end

@testitem "fitndgr for n = 2" begin
    using Statistics: mean
    x = sort(rand(10))
    y = @. 3 * x^2 + 2 * x + 1
    f = fitndgr(x, y, 2)
    @test f.R2 ≈ 1
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred == f.(x))

    x = Float32.(x)
    y = Float32.(y)
    f = fitndgr(x, y, 2)
    @test typeof(f.R2) == Float32
end

@testitem "fitndgr for n = 3" begin
    using Statistics: mean
    x = sort(rand(10))
    y = @. 4 * x^3 + 3 * x^2 + 2 * x + 1
    f = fitndgr(x, y, 3)
    @test f.R2 ≈ 1
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred == f.(x))

    x = Float32.(x)
    y = Float32.(y)
    f = fitndgr(x, y, 3)
    @test typeof(f.R2) == Float32
end

@testitem "fitndgr for n = 4" begin
    using Statistics: mean
    x = sort(rand(10))
    y = @. 6 * x^4 + 4 * x^3 + 3 * x^2 + 2 * x + 1
    f = fitndgr(x, y, 4)
    @test f.R2 ≈ 1
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred == f.(x))

    f = fitndgr(x, y, 4, l=[-Inf, -Inf, -Inf, 5.0, -Inf])
    @test f.R2 ≈ 1.0 atol = 1e-3
    @test f.lscoeff[4] ≈ 5.0 atol = 1e-5

    f = fitndgr(x, y, 4, u=[+Inf, +Inf, +Inf, -5.0, +Inf])
    @test f.R2 ≈ 1.0 atol = 1e-3
    @test f.lscoeff[4] ≈ -5.0 atol = 1e-5

    x = Float32.(x)
    y = Float32.(y)
    f = fitndgr(x, y, 4)
    @test typeof(f.R2) == Float32
end
