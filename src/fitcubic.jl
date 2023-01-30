#
# Cubic fit
#
struct Cubic{T} <: Fit{T}
    a::T
    b::T
    c::T
    d::T
    R::T
    x::Vector{T}
    y::Vector{T}
    ypred::Vector{T}
    residues::Vector{T}
end

"""
    fitcubic(x,y)

Obtains the cubic polynomial fit: ``y = ax^3 + bx^2 + cx + d``

Optional lower and upper bounds for a, b, and c can be provided using, for example:

```
fitcubic(x,y, l=lower(b=0.), u=upper(a=5.))
```

and `d` can be set to constant with, for example:

```
fitcubic(x,y,d=5.)
```

# Examples
```jldoctest
julia>  x = sort(rand(10)); y = x.^3 .+ rand(10);

julia> fit = fitcubic(x,y)

 ------------------- Cubic Fit ----------------- 

 Equation: y = ax^3 + bx^2 + cx + d 

 With: a = 12.637633791600711
       b = -19.648194970330454
       c = 10.018385827387148
       d = -0.8740912356800155

 Pearson correlation coefficient, R = 0.7831345513024988
 Average square residue = 0.0781543071776559

 Predicted Y: ypred = [0.24999805379642903, 0.3001612840610868...
 residues = [0.2238223147726266, 0.12656861200050698...

 ----------------------------------------------- 

```
"""
function fitcubic(
    X::AbstractArray{T}, Y::AbstractArray{T};
    l::lower=lower(), u::upper=upper(), d=nothing, options::Options=Options()
) where {T<:Real}
    # Check data
    X, Y = checkdata(X, Y, options)
    # Set bounds
    vars = [VarType(:a, Number, 1),
        VarType(:b, Number, 1),
        VarType(:c, Number, 1),
        VarType(:d, Nothing, 1)]
    lower, upper = setbounds(vars, l, u, T)
    if isnothing(d)
        # Set model
        @. model(x, p) = p[1] * x^3 + p[2] * x^2 + p[3] * x + p[4]
        # Fit
        fit = find_best_fit(model, X, Y, length(vars), options, lower, upper)
        # Analyze results and return
        R = pearson(X, Y, model, fit)
        x, y, ypred = finexy(X, options.fine, model, fit)
        return Cubic(fit.param..., R, x, y, ypred, fit.resid)
    else
        lower = lower[1:length(vars)-1]
        upper = upper[1:length(vars)-1]
        # Set model
        @. model_const(x, p) = p[1] * x^3 + p[2] * x^2 + p[3] * x + d
        # Fit
        fit = find_best_fit(model_const, X, Y, length(vars) - 1, options, lower, upper)
        # Analyze results and return
        R = pearson(X, Y, model_const, fit)
        x, y, ypred = finexy(X, options.fine, model, fit)
        return Cubic(fit.param..., d, R, x, y, ypred, fit.resid)
    end
end

"""
    (fit::EasyFit.Cubic)(x::Real)

Calling the the fitted estimator on a new sample generates point predictions. To compute predictions for multiple new data points, use broadcasting.

# Examples

```jldoctest
julia> x = sort(rand(10)); y = x.^3 .+ rand(10);

julia> f = fitcubic(x,y)

 ------------------- Cubic Fit ----------------- 

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
function (fit::EasyFit.Cubic)(x::Real)
    a = fit.a
    b = fit.b
    c = fit.c
    d = fit.d
    return a * x^3 + b * x^2 + c * x + d
end

function Base.show(io::IO, fit::Cubic)
    println(io,
        """
        ------------------- Cubic Fit -----------------

        Equation: y = ax^3 + bx^2 + cx + d

        With: a = $(fit.a)
              b = $(fit.b)
              c = $(fit.c)
              d = $(fit.d)

        Pearson correlation coefficient, R = $(fit.R))
        Average square residue = $(mean(fit.residues .^ 2))

        Predicted Y: ypred = [ $(fit.ypred[1]), $(fit.ypred[2]), ...]
        residues = [ $(fit.residues[1]), $(fit.residues[2]), ...]

        -----------------------------------------------"""
    )
end

export fitcubic

@testitem "fitcubic" begin
    using Statistics: mean
    x = sort(rand(10))
    y = @. 4 * x^3 + 3 * x^2 + 2 * x + 1
    f = fitcubic(x, y)
    @test f.R ≈ 1
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R^2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred == f.(x))

    x = Float32.(x)
    y = Float32.(y)
    f = fitcubic(x, y)
    @test typeof(f.R) == Float32
end
