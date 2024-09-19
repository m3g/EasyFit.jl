#
# Quadratic fit
# 

@kwdef struct Quadratic{TA,TB,TR,TX,TY}
    a::TA
    b::TB
    c::TY
    R::TR
    x::Vector{TX}
    y::Vector{TY}
    ypred::Vector{TY}
    residues::Vector{TY}
end

"""
    fitquad(x,y)
    fitquadratic(x,y)

Obtains the quadratic fit: ``y = a*x^2 + b*x + c``

Optional lower and upper bounds for a and b can be provided using, for example:

```
fitquad(x,y, lower(b=0.), upper(a=5.,b=7.) )
```
and the intercept `c` can be fixed with

```
fitquad(x,y, c=3.)
```

# Examples
```jldoctest
julia>  x = sort(rand(10)); y = x.^2 .+ rand(10);

julia> fit = fitquad(x,y)

 ------------------- Quadratic Fit ------------- 

 Equation: y = ax^2 + bx + c 

 With: a = 1.9829681649993036
       b = -1.24215737650827
       c = 0.9410816080128867

 Pearson correlation coefficient, R = 0.8452759310204063
 Average square residue = 0.039620067833833005

 Predicted Y: ypred = [0.778952191090992, 0.7759243614999851...
 residues = [0.0550252612868799, -0.15207394277809727...

 ----------------------------------------------- 

```
"""
function fitquadratic(
    X::AbstractArray{TX}, Y::AbstractArray{TY};
    l::lower=lower(), u::upper=upper(), c=nothing,
    options::Options=Options()
) where {TX<:Number,TY<:Number}
    onex = oneunit(TX)
    oney = oneunit(TY)
    # Check data
    X, Y, data_type = checkdata(X, Y, options)
    # Set bounds
    vars = [VarType(:a, Number, 1), VarType(:b, Number, 1), VarType(:c, Nothing, 1)]
    lower, upper = setbounds(vars, l, u, oney)
    if isnothing(c)
        @. model(x, p) = p[1] * x^2 + p[2] * x + p[3]
        fit = find_best_fit(model, X, Y, length(vars), options, lower, upper)
        R = pearson(X, Y, model, fit)
        x, y, ypred = finexy(X, options.fine, model, fit)
    else
        if unit(c) != unit(TY)
            throw(ArgumentError("The interecept c must have the same units as the dependent variable y: $(unit(TY))"))
        end
        c = ustrip(c)
        lower = lower[1:length(vars)-1]
        upper = upper[1:length(vars)-1]
        @. model_const(x, p) = p[1] * x^2 + p[2] * x + c
        fit = find_best_fit(model_const, X, Y, length(vars) - 1, options, lower, upper)
        R = pearson(X, Y, model_const, fit)
        x, y, ypred = finexy(X, options.fine, model_const, fit)
    end
    return Quadratic(
        a=fit.param[1] * oney / (onex^2),
        b=fit.param[2] * oney / onex,
        c=isnothing(c) ? fit.param[3] * oney : oney * c,
        R=R * onex * oney,
        x=x * onex,
        y=y * oney,
        ypred=ypred * oney,
        residues=fit.resid * oney
    )
end
const fitquad = fitquadratic

"""
    (fit::Quadratic)(x::Real)

Calling the the fitted estimator on a new sample generates point predictions. To compute predictions for multiple new data points, use broadcasting.

# Examples

```jldoctest
julia> x = sort(rand(10)); y = x.^2 .+ rand(10);

julia> f = fitquad(x,y)

 ------------------- Quadratic Fit ------------- 

 Equation: y = ax^2 + bx + c 

 With: a = 3.0661527272135043
       b = -1.2832262361743607
       c = 0.3650565332863989

 Pearson correlation coefficient, R = 0.8641384358642901
 Average square residue = 0.06365720683818799

 Predicted Y: ypred = [0.2860912283436672, 0.2607175680542409...
 residues = [0.11956927540814413, -0.26398094542690925...

 ----------------------------------------------- 


julia> f.(rand(10))
10-element Vector{Float64}:
 1.7769193832294496
 0.2489307162532048
 0.33127545070267345
 0.4422093927705098
 0.2974933484569105
 0.6299254836558978
 0.24575582331233475
 0.9185494116516484
 1.4615291776107249
 1.600204246377446
```
"""
function (fit::Quadratic)(x::Real)
    a = fit.a
    b = fit.b
    c = fit.c
    return a * x^2 + b * x + c
end

function Base.show(io::IO, fit::Quadratic)
    println(io, chomp("""
        ------------------- Quadratic Fit -------------

        Equation: y = ax^2 + bx + c

        With: a = $(fit.a)
              b = $(fit.b)
              c = $(fit.c)

        Pearson correlation coefficient, R = $(fit.R)
        Average square residue = $(mean(fit.residues .^ 2))

        Predicted Y: ypred = [$(fit.ypred[1]), $(fit.ypred[2]), ...]
        residues = [$(fit.residues[1]), $(fit.residues[2]), ...]

        -----------------------------------------------"""
   ))
end

export fitquad, fitquadratic

@testitem "fitquadratic" begin
    using Unitful
    using Statistics: mean
    x = sort(rand(10))
    y = @. 3x^2 + 2x + 1
    f = fitquadratic(x, y)
    @test f.R ≈ 1
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R^2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred == f.(x))

    x = Float32.(x)
    y = Float32.(y)
    f = fitquadratic(x, y)
    @test typeof(f.R) == Float32

    # with units
    x = sort(rand(10))u"s"
    xu = ustrip(x)
    y = (@. 3xu^2 + 2xu + 1)u"m"

    f = fitquadratic(x, y)
    @test f.R ≈ 1u"m*s"
    @test f.a ≈ 3u"m/s^2"
    @test f.b ≈ 2u"m/s"
    @test f.c ≈ 1u"m"

    f = fitquadratic(x, y; c=1u"m")
    @test f.R ≈ 1u"m*s"
    @test f.a ≈ 3u"m/s^2"
    @test f.b ≈ 2u"m/s"
    @test f.c ≈ 1u"m"
    @test_throws ArgumentError fitquadratic(x, y; c=1u"s")

    f = fitquadratic(x, y; l=lower(a=4.0))
    @test f.a ≈ 4u"m/s^2"
    f = fitquadratic(x, y; u=upper(a=0.0))
    @test f.a ≈ 0u"m/s^2"
    f = fitquadratic(x, y; l=lower(b=2.0))
    @test f.b ≈ 2u"m/s"
    f = fitquadratic(x, y; u=upper(b=-1.0))
    @test f.b ≈ -1u"m/s"
    f = fitquadratic(x, y; l=lower(a=5.0, b=3.0))
    @test f.a ≈ 5u"m/s^2"
    @test f.b ≈ 3u"m/s"
    f = fitquadratic(x, y; u=upper(a=-2.0, b=-1.0))
    @test f.a ≈ -2u"m/s^2"
    @test f.b ≈ -1u"m/s"
end
