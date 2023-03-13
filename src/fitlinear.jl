# 
# Linear fit
#

struct Linear{T} <: Fit{T}
    a::T
    b::T
    R::T
    x::Vector{T}
    y::Vector{T}
    ypred::Vector{T}
    residues::Vector{T}
end

"""
    fitlinear(x,y)

Obtains the linear fit: ``y = a*x + b``

Optional lower and upper bounds for a, and constant b can be provided using, for example:

```
fitlinear(x,y, l=lower(a=0.), b=3)
```

# Examples
```jldoctest
julia> x = sort(rand(10)) ; y = sort(rand(10));

julia> fit = fitlinear(x,y)

------------------- Linear Fit ------------- 

Equation: y = ax + b 

With: a = 1.0448783208110997
      b = 0.18817627115683894

Pearson correlation coefficient, R = 0.8818586822210751
Average absolute residue = 0.14274752107157443

Predicted Y: ypred = [0.1987357699444139, 0.32264343301109627...
residues = [0.1613987313816987, 0.22309410865095275...

-------------------------------------------- 

```
"""
function fitlinear(
    X::AbstractArray{T1}, Y::AbstractArray{T2};
    l::lower=lower(), u::upper=upper(), b=nothing,
    options::Options=Options()
) where {T1<:Real, T2<:Real}
    # Check data
    X, Y, data_type = checkdata(X, Y, options)
    # Set bounds
    vars = [VarType(:a, Number, 1), VarType(:b, Nothing, 1)]
    lower, upper = setbounds(vars, l, u, data_type)
    if isnothing(b)
        # Set model
        @. model(x, p) = p[1] * x + p[2]
        # Initial point
        p0 = Vector{data_type}(undef, 2)
        initP!(p0, options, lower, upper)
        # Fit
        fit = curve_fit(model, X, Y, p0, lower=lower, upper=upper)
        # Analyze results and return
        R = pearson(X, Y, model, fit)
        x, y, ypred = finexy(X, length(X), model, fit)
        return Linear(fit.param..., R, x, y, ypred, fit.resid)
    else
        lower = [lower[1]]
        upper = [upper[1]]
        # Set model
        @. model_const(x, p) = p[1] * x + b
        # Initial point
        p0 = Vector{data_type}(undef, 1)
        initP!(p0, options, lower, upper)
        # Fit
        fit = curve_fit(model_const, X, Y, p0, lower=lower, upper=upper)
        # Analyze results and return
        R = pearson(X, Y, model_const, fit)
        x, y, ypred = finexy(X, length(X), model_const, fit)
        return Linear(fit.param..., b, R, x, y, ypred, fit.resid)
    end
end

"""
    (fit::Linear)(x::Real)

Calling the the fitted estimator on a new sample generates point predictions. To compute predictions for multiple new data points, use broadcasting.

# Examples

```jldoctest
julia> x = sort(rand(10)) ; y = sort(rand(10));

julia> f = fitlinear(x,y)

 ------------------- Linear Fit ------------- 

 Equation: y = ax + b 

 With: a = 0.7492056732121763
       b = 0.19051493263850805

 Pearson correlation coefficient, R = 0.9880617647915537
 Average square residue = 0.0011903365407044974

 Predicted Y: ypred = [0.19131831893286483, 0.2588265305624418...
 residues = [-0.015828020422002875, -0.0503384398812427...

 -------------------------------------------- 


julia> f.(rand(10))
10-element Vector{Float64}:
 0.318112972601176
 0.5541324942065607
 0.2684448170646049
 0.2448341076856998
 0.19914167794590798
 0.7043365726554222
 0.47993606210138606
 0.6726328561177188
 0.3094063592157996
 0.40113908380656205
```
"""
function (fit::Linear)(x::Real)
    a = fit.a
    b = fit.b
    return a * x + b
end

function Base.show(io::IO, fit::Linear)
    println(io,
        """
        ------------------- Linear Fit -------------

        Equation: y = ax + b

        With: a = $(fit.a)
              b = $(fit.b)

        Pearson correlation coefficient, R = $(fit.R)
        Average square residue = $(mean(fit.residues .^ 2))

        Predicted Y: ypred = [$(fit.ypred[1]), $(fit.ypred[2]), ...]
        residues = [$(fit.residues[1]), $(fit.residues[2]), ...]

        --------------------------------------------"""
    )
end

export fitlinear

@testitem "fitlinear" begin
    using Statistics: mean
    x = sort(rand(10))
    y = @. 2x + 1
    f = fitlinear(x, y)
    @test f.R ≈ 1
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R^2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred == f.(x))

    x = Float32.(x)
    y = Float32.(y)
    f = fitlinear(x, y)
    @test typeof(f.R) == Float32
end
