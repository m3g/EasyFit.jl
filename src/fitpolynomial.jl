#
# Polynomial fit
#

struct Polynomial{T} <: Fit{T}
    θ::Vector{T}
    d::T
    R::T
    x::Vector{T}
    y::Vector{T}
    ypred::Vector{T}
    residues::Vector{T}
end

"""
    fitpolynomial(x,y,n)

Obtains the `n`-polynomial fit: ``y = Σᵢ θ[i] xⁱ + d``

Optional lower and upper bounds for a, b, and c can be provided using, for example:

```
fitpolynomial(x,y, l=lower(b=0.), u=upper(a=5.))
```

and `d` can be set to constant with, for example:

```
fitpolynomial(x,y,n,d=5.)
```

# Examples
```jldoctest
julia>  x = sort(rand(10)); y = x.^3 .+ rand(10);

```
"""
function fitpolynomial(X::AbstractArray{<:Real}, Y::AbstractArray{<:Real}, n::Int=1;
    l::lower=lower(), u::upper=upper(), d=nothing,
    options::Options=Options())
    # Check data
    X, Y = checkdata(X, Y, options)
    # Set bounds
    vars = [VarType(:a, Vector, n), VarType(:d, Nothing, 1)]
    lower, upper = setbounds(vars, l, u)
    if isnothing(d)
        # Set model
        @. model(x, p) = permutedims([x^i for i in 1:n]) * p
        # Fit
        fit = find_best_fit(model, X, Y, length(vars), options, lower, upper)
        # Analyze results and return
        R = pearson(X, Y, model, fit)
        x, y, ypred = finexy(X, options.fine, model, fit)
        return Polynomial(fit.param..., R, x, y, ypred, fit.resid)
    else
        lower = lower[1:length(vars)-1]
        upper = upper[1:length(vars)-1]
        # Set model
        @. model_const(x, p) = permutedims([x^i for i in 1:n]) * p + d
        # Fit
        fit = find_best_fit(model_const, X, Y, length(vars) - 1, options, lower, upper)
        # Analyze results and return
        R = pearson(X, Y, model_const, fit)
        x, y, ypred = finexy(X, options.fine, model, fit)
        return Polynomial(fit.param..., d, R, x, y, ypred, fit.resid)
    end
end

"""
    (fit::EasyFit.Polynomial)(x::Real)

Calling the the fitted estimator on a new sample generates point predictions. To compute predictions for multiple new data points, use broadcasting.

# Examples

```jldoctest
julia> x = sort(rand(10)); y = x.^3 .+ rand(10);

julia> f = fitpolynomial(x,y)

 ------------------- Polynomial Fit ----------------- 

 Equation: y = ax^3 + bx^2 + cx + d 

 With: a = 3.498571133673037
       b = -5.75292789995513
       c = 2.626129810011887
       d = 0.6361773562878126

 Pearson correlation coefficient, R = 0.7405690253097572
 Average square residue = 0.01215483592609077

 Predicted Y: ypred = [0.6416314330095221, 0.6417874373639705...
 residues = [-0.13182717628179608, -0.01592993507117535...

 ----------------------------------------------- 


julia> f.(rand(10))
10-element Vector{Float64}:
 0.8761239348448231
 0.9115358893542463
 0.9121562305431836
 0.8919530945018805
 ⋮
 0.81693749334824
 0.9622975666245418
 0.9753695182250022
```
"""
function (fit::EasyFit.Polynomial)(x::Real)
    θ = fit.θ
    d = fit.d
    return permutedims([x^i for i in 1:n]) * θ + d
end

function Base.show(io::IO, fit::Polynomial)
    println(io,"""
    ------------------- Polynomial Fit -----------------

    Equation: y = Σᵢ θ[i] xⁱ + d

    Pearson correlation coefficient, R = $(fit.R)
    Average square residue = $(mean(fit.residues .^ 2))

    Predicted Y: ypred = [$(fit.ypred[1]), $(fit.ypred[2]), ...]
    residues = [$(fit.residues[1]), $(fit.residues[2]), ...]

    -----------------------------------------------
    """)
end

export fitpolynomial

@testitem "fitpolynomial" begin
    using Statistics: mean
    x = sort(rand(10))
    y = @. 3x^2 + 2x + 1
    f = fitpolynomial(x, y, 3)
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R^2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred == f.(x))
end
