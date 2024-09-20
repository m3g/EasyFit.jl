# 
# Linear fit
#

@kwdef struct Linear{TA,TR,TX,TY}
    a::TA
    b::TY
    R::TR
    x::Vector{TX}
    y::Vector{TY}
    ypred::Vector{TY}
    residues::Vector{TY}
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
        R = pearson(X, Y, model, fit)
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
        R = pearson(X, Y, model_const, fit)
        x, y, ypred = finexy(X, length(X), model_const, fit)
    end
    return Linear(
        a=fit.param[1] * oney / onex,
        b=isnothing(b) ? fit.param[2] * oney : oney * b,
        R=R .* onex * oney,
        x=x .* onex,
        y=y .* oney,
        ypred=ypred .* oney,
        residues=fit.resid .* oney
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

        With: a = $(fit.a)
              b = $(fit.b)

        Pearson correlation coefficient, R = $(fit.R)
        Average square residue = $(mean(fit.residues .^ 2))

        Predicted Y: ypred = [$(fit.ypred[1]), $(fit.ypred[2]), ...]
        residues = [$(fit.residues[1]), $(fit.residues[2]), ...]

        --------------------------------------------
        """
    ))
end

export fitlinear

@testitem "fitlinear" begin
    using Statistics: mean
    using Unitful
    x = sort(rand(10))
    y = @. 2x + 1
    f = fitlinear(x, y)
    @test f.R ≈ 1
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R^2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred == f.(x))

    # with units
    x = sort(rand(10))u"s"
    y = (@. 2(ustrip(x)) + 1)u"m"
    f = fitlinear(x, y)
    @test f.R ≈ 1u"m*s"
    @test f.a ≈ 2u"m/s"
    @test f.b ≈ 1u"m"
    @test f(1.0u"s") ≈ 3.0u"m"

    f = fitlinear(x, y; b=1u"m")
    @test f.R ≈ 1u"m*s"
    @test f.a ≈ 2u"m/s"
    @test f.b ≈ 1u"m"
    @test_throws ArgumentError fitlinear(x, y; b=1u"s")

    f = fitlinear(x, y; l=lower(a=0.0))
    @test f.R ≈ 1u"m*s"
    @test f.a ≈ 2u"m/s"
    @test f.b ≈ 1u"m"

    f = fitlinear(x, y; u=upper(a=3.0))
    @test f.R ≈ 1u"m*s"
    @test f.a ≈ 2u"m/s"
    @test f.b ≈ 1u"m"
end

