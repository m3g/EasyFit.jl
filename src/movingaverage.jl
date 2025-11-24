#
# Computes the moving average of the data
#

struct MovingAverage{T} <: Fit{T}
    n::Int
    x::Vector{T}
    R2::T
    residues::Vector{T}
end

"""
    movingaverage(x,n)
    movavg(x,n)

Computes a moving average of `x[i]` in the range `i±(n-1)/2`. If `n` is even, we do `n ←  n + 1`.

# Example

```jldoctest
julia> x = rand(10);

julia> movingaverage(x,3)

 ------------------- Moving Average ----------

 Number of points averaged: 3 (± 1 points)

 Correlation coefficient, R² = 0.3532754137104625

 Averaged X: x = [0.5807828672543551, 0.40496733381946143...
 residues = [-0.22791917753944557, 0.4037347109743393...

 --------------------------------------------
```

"""
function movingaverage(X::AbstractArray{T}, n::Int) where {T<:Real}
    if n == 0
        error(" Please set n with movavg(x,n=10), for example. ")
    end
    if !isodd(n)
        n = n + 1
    end
    delta = div(n - 1, 2)
    nX = length(X)
    x = similar(X)
    for i in 1:nX
        jmin = max(i - delta, 1)
        jmax = min(i + delta, nX)
        x[i] = zero(T)
        for j in jmin:jmax
            x[i] = x[i] + X[j]
        end
        x[i] = x[i] / (jmax - jmin + 1)
    end
    residues = similar(X)
    for i in 1:nX
        residues[i] = X[i] - x[i]
    end
    R = Statistics.cor(x, X)
    return MovingAverage(n, x, R, residues)
end
movingaverage(X::AbstractArray{<:Real}; n::Int=0) = movingaverage(X, n)
const movavg = movingaverage

function Base.show(io::IO, fit::MovingAverage)
    println(io,
        """
        ------------------- Moving Average ----------

        Number of points averaged: $(fit.n) (± $(round(Int,(fit.n-1)/2)) points.

        Correlation coefficient, R² = $(fit.R2)

        Averaged X: x = [$(fit.x[1]), $(fit.x[2]), ...]
        residues = [$(fit.residues[1]), $(fit.residues[2]), ...]

        --------------------------------------------"""
    )
end

export movingaverage, movavg

@testitem "moving average" begin
    using ShowMethodTesting
    x = sin.(-5:0.13:5)
    fit = movingaverage(x,n=10)
    @test eltype(fit.x) == Float64
    @test parse_show(fit) ≈ """
    ------------------- Moving Average ----------
    
    Number of points averaged: 11 (± 5 points.
    
    Correlation coefficient, R² = 0.9997773413109331
    
    Averaged X: x = [0.9748470959093063, 0.9614698001834519, ...]
    residues = [-0.01592282124616784, 0.026135273736763498, ...]
    
    --------------------------------------------
    """ float_match = (x,y) -> isapprox(x,y; atol=1e-3)

    x = Float32.(x)
    fit = movingaverage(x,n=10)
    @test eltype(fit.x) == Float32
end

