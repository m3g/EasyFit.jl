#
# Exponential fit
#

struct SingleExponential{T} <: Fit{T}
    a::T
    b::T
    c::T
    R2::T
    x::Vector{T}
    y::Vector{T}
    ypred::Vector{T}
    residues::Vector{T}
end

struct MultipleExponential{T} <: Fit{T}
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

function sum_of_exps(x::Real, p::AbstractVector)
    n = div(length(p) - 1, 2)
    f = p[length(p)]
    for i in 1:n
        f = f + p[i] * exp(-x / p[n+i])
    end
    return f
end

function exp_model(X::AbstractVector{T}, p::AbstractVector{T}) where {T}
    f = Vector{T}(undef, length(X))
    for i in eachindex(X)
        f[i] = sum_of_exps(X[i], p)
    end
    return f
end

function sum_of_exps_const(x::Real, p::AbstractVector, c)
    n = div(length(p), 2)
    f = c
    for i in 1:n
        f = f + p[i] * exp(-x / p[n+i])
    end
    return f
end

function exp_model_const(X::AbstractArray{T}, p::AbstractVector{T}, c) where {T<:Real}
    f = Vector{T}(undef, length(X))
    for i in eachindex(X)
        f[i] = sum_of_exps_const(X[i], p, c)
    end
    return f
end

"""
    fitexponential(x,y; n::Int=1)
    fitexp(x,y; n::Int=1)
  
Obtains single or multiexponential fits: ``y = a*exp(-x/b) + c`` or  ``y = sum(a[i]*exp(-x/b[i]) for i in 1:N) + c``

Lower and upper bounds can be optionall set, and the intercept `c` can set to be constant:

For single exponentials all parameters are scalars:

```
fitexp(x,y,l=lower(b=0.),u=upper(b=10.),c=5.)
```

For multiple exponentials, `a` and `b` bounds must be vectors of dimension `N`.

```
fitexp(x,y,n=2,l=lower(a=[0.,0.]),u=upper(b=[-100.,-5.]))
```

# Examples
```jldoctest
julia> x = sort(rand(10)); y = rand()*exp.(sort(rand(10)));

julia> fit = fitexp(x,y,l=lower(b=[0.,0.]),n=2)

 -------- Multiple-exponential fit -------------

 Equation: y = sum(a[i] exp(-x/b[i]) for i in 1:2) + c

 With: a = [6.60693727987886e-13, 0.6249999999993409]
       b = [0.02688289803014393, 0.5000000000002596]
       c = 0.37499999999999856

 Correlation coefficient, R² = 1.0
 Average square residue = 1.1639900380979497e-29

 Predicted Y: ypred = [1.0000000000000002, 0.4595845520228801...
 residues = [2.220446049250313e-16, -2.831068712794149e-15...

 -----------------------------------------------

```
"""
function fitexponential(
    X::AbstractArray{TX}, Y::AbstractArray{TY};
    n::Int=1, c=nothing, l::lower=lower(), u::upper=upper(), options::Options=Options()
) where {TX<:Real,TY<:Real}
    onex = oneunit(TX)
    oney = oneunit(TY)
    # Check data
    X, Y, data_type = checkdata(X, Y, options)
    if isnothing(c)
        # Set model
        model(x, p) = exp_model(x, p)
        # Number of exponentials
        # Set bounds
        if n == 1
            vars = [VarType(:a, Number, 1), VarType(:b, Number, 1), VarType(:c, Nothing, 1)]
            lower, upper = setbounds(vars, l, u, oney)
        else
            vars = [VarType(:a, Vector, n), VarType(:b, Vector, n), VarType(:c, Nothing, 1)]
            lower, upper = setbounds(vars, l, u, oney)
        end
        # Fit
        fit = find_best_fit(model, X, Y, 2 * n + 1, options, lower, upper)
        # Analyze and return
        R = R2(X, Y, model, fit)
        x, y, ypred = finexy(X, options.fine, model, fit)
        if n == 1
            return SingleExponential(fit.param[1], fit.param[2], fit.param[3],
                R, x, y, ypred, fit.resid)
        else
            ind = collect(1:n)
            sort!(ind, by=i -> fit.param[n+i])
            a = fit.param[1:n][ind]
            b = fit.param[n+1:2*n][ind]
            return MultipleExponential(n, a, b, fit.param[2*n+1],
                R, x, y, ypred, fit.resid)
        end
    else
        # Set model
        model_const(x, p) = exp_model_const(x, p, c)
        # Number of exponentials
        # Set bounds
        if n == 1
            vars = [VarType(:a, Number, 1), VarType(:b, Number, 1), VarType(:c, Nothing, 1)]
            lower, upper = setbounds(vars, l, u, oney)
        else
            vars = [VarType(:a, Vector, n), VarType(:b, Vector, n), VarType(:c, Nothing, 1)]
            lower, upper = setbounds(vars, l, u, oney)
        end
        lower = lower[1:2*n]
        upper = upper[1:2*n]
        # Fit
        fit = find_best_fit(model_const, X, Y, 2 * n, options, lower, upper)
        # Analyze and return
        R = R2(X, Y, model_const, fit)
        x, y, ypred = finexy(X, options.fine, model_const, fit)
        if n == 1
            return SingleExponential(fit.param[1], fit.param[2], c, R, x, y, ypred, fit.resid)
        else
            ind = collect(1:n)
            sort!(ind, by=i -> fit.param[n+i])
            a = fit.param[1:n][ind]
            b = fit.param[n+1:2*n][ind]
            return MultipleExponential(n, a, b, c, R, x, y, ypred, fit.resid)
        end
    end
end
fitexp = fitexponential

"""
    (fit::SingleExponential)(x::Real)

Calling the the fitted estimator on a new sample generates point predictions. To compute predictions for multiple new data points, use broadcasting.

# Examples

```jldoctest
julia> x = sort(rand(10)); y = rand()*exp.(sort(rand(10)));

julia> f = fitexp(x,y,l=lower(b=0.),u=upper(b=10.),c=5.)

 ------------ Single Exponential fit ----------- 

 Equation: y = a exp(-x/b) + c

 With: a = -4.663077750813696
       b = 3.503891711388752
       c = 5.0

 Correlation coefficient, R² = 0.8451556303228667
 Average square residue = 0.018962968678623245

 Predicted Y: ypred = [0.5349351918795522, 0.8336327629821874...
 residues = [-0.1770144709939867, 0.026331423737706694...

 ----------------------------------------------- 


julia> f.(rand(10))
10-element Vector{Float64}:
 1.3311695732688578
 0.3472845242931859
 1.4590115873141651
 0.3763535792927408
 0.5756106398622487
 1.423007971271891
 0.8893556848381099
 0.7803804518815749
 0.9734992788718548
 0.5963345544654599
```
"""
function (fit::SingleExponential)(x::Real)
    a = fit.a
    b = fit.b
    c = fit.c
    return a * exp(-x / b) + c
end

"""

# Examples

```jldoctest
julia> x = sort(rand(10)); y = rand()*exp.(sort(rand(10)));

julia> f = fitexp(x,y,l=lower(b=[0.,0.]),n=2)

 -------- Multiple-exponential fit ------------- 

 Equation: y = sum(a[i] exp(-x/b[i]) for i in 1:2) + c 

 With: a = [70.52352486967449, -71.7586710966385]
       b = [0.4660083805047125, 0.4780084734705192]
       c = 1.8482858404652986

 Correlation coefficient, R² = 0.9933941403034768
 Average square residue = 0.0007163482924343316

 Predicted Y: ypred = [0.5864990290966858, 0.5582571702322241...
 residues = [0.003833387360068774, -0.026186725644245512...

 ----------------------------------------------- 


julia> f.(rand(10))
10-element Vector{Float64}:
 0.7484563114305345
 0.9716494335236276
 0.5700590449540968
 1.2075688838563046
 1.2280682982441256
 0.5361263670282497
 0.6277801581255005
 0.7234067850308292
 0.5698309060534725
 0.5441089815268014
```
"""
function (fit::MultipleExponential)(x::Real)
    a = fit.a
    b = fit.b
    c = fit.c
    n = length(a)
    return sum(a[i] * exp(-x / b[i]) for i in 1:n) + c
end

function Base.show(io::IO, fit::SingleExponential)
    println(io,
        """
        ------------ Single Exponential fit -----------

        Equation: y = a exp(-x/b) + c

        With: a = $(fit.a)
              b = $(fit.b)
              c = $(fit.c)

        Correlation coefficient, R² = $(fit.R2)
        Average square residue = $(mean(fit.residues .^ 2))

        Predicted Y: ypred = [ $(fit.ypred[1]), $(fit.ypred[2]), ...]
        residues = [ $(fit.residues[1]), $(fit.residues[2]), ...]

        -----------------------------------------------"""
    )
end

function Base.show(io::IO, fit::MultipleExponential)
    println(io,
        """
        -------- Multiple-exponential fit -------------

        Equation: y = sum(a[i] exp(-x/b[i]) for i in 1:$(fit.n)) + c

        With: a = $(fit.a)
              b = $(fit.b)
              c = $(fit.c)

        Correlation coefficient, R² = $(f.R2)
        Average square residue = $(mean(fit.residues .^ 2))

        Predicted Y: ypred = [$(fit.ypred[1]), $(fit.ypred[2]), ...]
        residues = [$(fit.residues[1]), $(fit.residues[2]), ...]

        -----------------------------------------------"""
    )
end

export fitexp, fitexponential

@testitem "fitexponential" begin
    using ShowMethodTesting
    using Statistics: mean
    x = 0.1:0.13:5
    y = @. 3 * exp(-x / 2) + 1
    f = fitexp(x, y)
    @test f.R2 > 0.9
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred ≈ f.(x))
    @test parse_show(f) ≈ """
    ------------ Single Exponential fit -----------

    Equation: y = a exp(-x/b) + c
    
    With: a = 3.0
          b = 2.0
          c = 1.0
    
    Correlation coefficient, R² = 1.0
    Average square residue = 3.8924057823405186e-33
    
    Predicted Y: ypred = [ 3.853688273502142, 3.674098431720494, ...]
    residues = [ 0.0, 0.0, ...]
    
    -----------------------------------------------
    """ float_match = (x,y) -> isapprox(x,y; atol=1e-3)

    f = fitexp(x, y; n=2)
    @test f.R2 > 0.9
    @test all(f.ypred - y .== f.residues)
    ss_res = sum(f.residues .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    @test isapprox(f.R2, 1 - (ss_res / ss_tot), atol=1e-5)
    @test all(f.ypred ≈ f.(x))

    x = Float32.(x)
    y = Float32.(y)
    f = fitexp(x, y)
    @test typeof(f.R2) == Float32

    f = fitexp(x, y; n=2)
    @test typeof(f.R2) == Float32
end
