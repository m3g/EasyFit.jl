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

Use `step` to control the bin step. Use `norm=(0 or 1)` to set if the number of data points or the probability of finding a data point
within ± step/2 will be output.

By default, `norm=1` (probability) and the step is `(xmax-xmin)/100`.

# Arguments

- `x::AbstractVector{T}`: The data points to compute the density function from.
- `nbins::Int`: Number of bins to use for the density function.
- `step::T`: Step size for the bins. If not provided, it is calculated
    as `(xmax - xmin) / nbins`.
- `steptype::String`: Type of step size, either "absolute" or "relative". Default is "absolute".
- `vmin::T`: Minimum value for the bins. If not provided, it is calculated as the minimum of `x`.
- `vmax::T`: Maximum value for the bins. If not provided, it is calculated as the maximum of `x`.
- `norm::Int64`: Normalization type, either `0` (number of data points) or `1` (probability). Default is `1`.

# Returns

A `Density` object containing:

- `x`: Bin centers of the density function.
- `d`: Density values corresponding to the bin centers.
- `step`: The step size used for the bins.
- `norm`: Normalization type used (0 for number of data points, 1 for probability).

# Example
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

    x = zeros(T, nbins)
    df = zeros(T, nbins)
    binstep = (vmax - vmin) / nbins
    for i in 1:nbins
        x[i] = vmin + (i - 1) * binstep + binstep / 2
        nv = 0
        for j in 1:ndata
            if (v[j] > x[i] - step / 2) && (v[j] <= x[i] + step / 2)
                nv = nv + 1
            end
        end
        if norm == 0
            df[i] = nv
        elseif norm == 1
            df[i] = nv / ndata
        end
    end
    return Density(x, df, T(step), norm)
end
export fitdensity

function Base.show(io::IO, d::Density)
    if d.norm == 0
        s = "number of"
    end
    if d.norm == 1
        s = "probability of finding"
    end
    print(io,
        """
        ------------------- Density -------------

         Fields:
            `x` contains the bin centers of the density function
            `d` contains the $s data points within x ± $(round(d.step/2,sigdigits=3))

        -----------------------------------------"""
    )
end

@testitem "fitdensity" begin
    using ShowMethodTesting
    using Statistics: std

    x = 0:0.1:20
    fit = fitdensity(x)
    @test eltype(fit.d) == Float64
    @test parse_show(fit) ≈ """
    ------------------- Density -------------
    
     Fields:
        `x` contains the bin centers of the density function
        `d` contains the probability of finding data points within x ± 0.1
    
    -----------------------------------------
    """

    x = Float32.(x)
    fit = fitdensity(x)
    @test eltype(fit.d) == Float32
end
