#
# Computes the density function of a list of values
#
struct Density{T} <: Fit{T}
    x::Vector{T}
    d::Vector{T}
    step::T
    norm::Int
end

"""
    fitdensity(x; step, norm)

Obtains the density function given data sampled.

Use `step=(Float64)` to control the bin step. Use `norm=(0 or 1)` to set
if the number of data points or the probability of finding a data point
within ± step/2 will be output.

By default, `norm=1` (probability) and the step is `(xmax-xmin)/100`.

# Examples
```jldoctest
julia> x = randn(1000)

julia> d = fitdensity(x)

 ------------------- Density -------------

  `d` contains the number of data points within x ± 0.028325979340904375

 -----------------------------------------

```
"""
function fitdensity(
    v::AbstractVector{T};
    nbins=nothing,
    step=nothing,
    steptype="absolute",
    vmin=nothing,
    vmax=nothing,
    norm::Int64=1
) where {T<:Real}

    ndata = length(v)
    isnothing(vmin) && (vmin = minimum(v))
    isnothing(vmax) && (vmax = maximum(v))
    isnothing(nbins) && (nbins = 100)
    if isnothing(step)
        step = (vmax - vmin) / nbins
    else
        # By default, the step size is absolute
        if steptype == "relative"
            step = step * (vmax - vmin) / nbins
        elseif steptype != "absolute"
            throw(ArgumentError(" steptype must be \"relative\" or \"absolute\""))
        end
    end

    x = Vector{T}(undef, nbins)
    df = Vector{T}(undef, nbins)

    binstep = (vmax - vmin) / nbins
    for i in 1:nbins
        x[i] = vmin + (i - 1) * binstep + binstep / 2
        nv = 0
        for j in 1:ndata
            if (v[j] > x[i] - step / 2) && (v[j] <= x[i] + step / 2)
                nv = nv + 1
            end
        end
        binsize = min(vmax, x[i] + step / 2) - max(vmin, x[i] - step / 2)
        if norm == 0
            df[i] = nv
        elseif norm == 1
            df[i] = nv / (binsize * ndata)
        end
    end

    return Density(x, df, step, norm)
end
export fitdensity

function Base.show(io::IO, d::Density)
    if d.norm == 0
        s = "number of"
    end
    if d.norm == 1
        s = "probability of finding"
    end
    println(io,
        """
        ------------------- Density -------------

         d contains the $s data points within x ± $(round(d.step/2,sigdigits=3))

        -----------------------------------------"""
    )
end

@testitem "fitdensity" begin
    using Statistics: std

    x = randn(100)
    fit = fitdensity(x)
    @test eltype(fit.d) == Float64

    x = Float32.(x)
    fit = fitdensity(x)
    @test eltype(fit.d) == Float32
end
