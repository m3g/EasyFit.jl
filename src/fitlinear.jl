# 
# Linear fit
#

@kwdef struct Linear{TA,TR,TX,TY}
    a::TA
    b::TY
    R2::TR
    x::Vector{TX}
    y::Vector{TY}
    ypred::Vector{TY}
    residues::Vector{TY}
    sd_a::TA
    sd_b::TY
end

"""
    fitlinear(x,y)

Obtains the linear fit: ``y = a*x + b``

Optional lower and upper bounds for a, and constant b can be provided using, for example:

The output is a `Linear` struct with fields:
- `a`: slope of the linear fit
- `b`: intercept of the linear fit
- `R2`: R2 correlation coefficient
- `x`: x data used in the fit
- `y`: y data used in the fit
- `ypred`: predicted y values from the fit
- `residues`: residues of the fit (y - ypred)
- `sd_a`: standard deviation (error) in a
- `sd_b`: standard deviation (error) in b

```
fitlinear(x,y, l=lower(a=0.), b=3)
```

# Examples
```jldoctest
julia> x = sort(rand(10)) ; y = sort(rand(10));

julia> fit = fitlinear(x,y)
------------------- Linear Fit -------------

Equation: y = ax + b

With: a = 1.158930569179642 ± 0.6538824813074927
      b = -0.1251714526967127 ± 0.3588742142656203

Correlation coefficient, R² = 0.9696101474036224
Average square residue = 0.004113279571449428

Predicted Y: ypred = [0.1044876257374221, 0.2072397615587609, ...]
residues = [0.08428396483020295, 0.05555828441380739, ...]

--------------------------------------------

```
"""
function fitlinear(
    X::AbstractArray{TX}, Y::AbstractArray{TY};
    l::lower=lower(), u::upper=upper(), b=nothing,
    options::Options=Options()
) where {TX<:Number,TY<:Number}
    # Check units
    # Check data
    onex = oneunit(TX)
    oney = oneunit(TY)
    X, Y, data_type = checkdata(X, Y, options)
    # Set bounds
    vars = [VarType(:a, Number, 1), VarType(:b, Nothing, 1)]
    lower, upper = setbounds(vars, l, u, oney)
    if isnothing(b)
        @. model(x, p) = p[1] * x + p[2]
        p0 = zeros(data_type, 2)
        initP!(p0, options, lower, upper)
        fit = curve_fit(model, X, Y, p0, lower=lower, upper=upper)
        _R2 = R2(X, Y, model, fit)
        x, y, ypred = finexy(X, length(X), model, fit)
    else
        if unit(b) != unit(TY)
            throw(ArgumentError("The intercept b must have the same units as the dependent variable y: $(unit(TY))"))
        end
        b = ustrip(b)
        lower = [first(lower)]
        upper = [first(upper)]
        @. model_const(x, p) = p[1] * x + b
        p0 = zeros(data_type, 1)
        initP!(p0, options, lower, upper)
        fit = curve_fit(model_const, X, Y, p0, lower=lower, upper=upper)
        _R2 = R2(X, Y, model_const, fit)
        x, y, ypred = finexy(X, length(X), model_const, fit)
    end
    a=fit.param[1]
    b=isnothing(b) ? fit.param[2] : b
    _X_mean = mean(X)
    _X_var = sum((X[i] - _X_mean)^2 for i in eachindex(X))
    sd_a = sqrt(sum((Y[i] - (a*X[i] + b))^2 for i in eachindex(X)) / (_X_var * (length(X) - 2)))
    sd_b = sd_a * sqrt(mean(abs2, X))
    return Linear(
        a=a * oney / onex,
        b=b * oney,
        R2=_R2 .* onex * oney,
        x=x .* onex,
        y=y .* oney,
        ypred=ypred .* oney,
        residues=fit.resid .* oney,
        sd_a=sd_a * oney / onex,
        sd_b=sd_b * oney,
    )
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

With: a = 1.158930569179642 ± 0.6538824813074927
      b = -0.1251714526967127 ± 0.3588742142656203

Correlation coefficient, R² = 0.9696101474036224
Average square residue = 0.004113279571449428

Predicted Y: ypred = [0.1044876257374221, 0.2072397615587609, ...]
residues = [0.08428396483020295, 0.05555828441380739, ...]

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
function (fit::Linear)(x::Number)
    a = fit.a
    b = fit.b
    return a * x + b
end

function Base.show(io::IO, fit::Linear)
    println(io, chomp(
        """
        ------------------- Linear Fit -------------

        Equation: y = ax + b

        With: a = $(fit.a) ± $(fit.sd_a)
              b = $(fit.b) ± $(fit.sd_b)

        Correlation coefficient, R² = $(fit.R2)
        Average square residue = $(mean(fit.residues .^ 2))

        Predicted Y: ypred = [$(fit.ypred[1]), $(fit.ypred[2]), ...]
        residues = [$(fit.residues[1]), $(fit.residues[2]), ...]

        --------------------------------------------
        """
    ))
end

export fitlinear

@testitem "fitlinear" begin
    using ShowMethodTesting
    using Statistics: mean
    using Unitful
    x = sort(rand(10))
    y = @. 2x + 1
    f = fitlinear(x, y)
    @test f.R2 ≈ 1
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred == f.(x))

    # with units
    x = sort(rand(10))u"s"
    y = (@. 2(ustrip(x)) + 1)u"m"
    f = fitlinear(x, y)
    @test f.R2 ≈ 1u"m*s"
    @test f.a ≈ 2u"m/s"
    @test f.b ≈ 1u"m"
    @test f(1.0u"s") ≈ 3.0u"m"

    f = fitlinear(x, y; b=1u"m")
    @test f.R2 ≈ 1u"m*s"
    @test f.a ≈ 2u"m/s"
    @test f.b ≈ 1u"m"
    @test_throws ArgumentError fitlinear(x, y; b=1u"s")

    f = fitlinear(x, y; l=lower(a=0.0))
    @test f.R2 ≈ 1u"m*s"
    @test f.a ≈ 2u"m/s"
    @test f.b ≈ 1u"m"

    f = fitlinear(x, y; u=upper(a=3.0))
    @test f.R2 ≈ 1u"m*s"
    @test f.a ≈ 2u"m/s"
    @test f.b ≈ 1u"m"

    data = [
        1.47	52.21
        1.50	53.12
        1.52	54.48
        1.55	55.84
        1.57	57.20
        1.60	58.57
        1.63	59.93
        1.65	61.29
        1.68	63.11
        1.70	64.47
        1.73	66.28
        1.75	68.10
        1.78	69.92
        1.80	72.19
        1.83	74.46
    ]

    fit = fitlinear(data[:,1], data[:,2])
    @test fit.a ≈ 61.272 atol=1e-3
    @test fit.b ≈ -39.062 atol=1e-3
    @test fit.sd_a ≈ 1.776 atol=1e-3
    @test fit.sd_b ≈ 2.938 atol=1e-3
    @test fit.R2 ≈ 0.989 atol=1e-3

    @test parse_show(fit) ≈ """
    ------------------- Linear Fit -------------
    
    Equation: y = ax + b
    
    With: a = 61.27218654191589 ± 1.7759227522153571
          b = -39.0619559185216 ± 2.93800106718343
    
    Correlation coefficient, R² = 0.989196922445797
    Average square residue = 0.49937056025883503
    
    Predicted Y: ypred = [51.00815829809476, 52.846323894352246, ...]
    residues = [-1.2018417019052379, -0.2736761056477519, ...]
    
    --------------------------------------------
    """

end

