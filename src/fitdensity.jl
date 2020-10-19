#
# Computes the density function of a list of values
#
struct Density
  x :: Vector{Float64}
  d :: Vector{Float64}
  step :: Float64
  norm :: Int64
end

"""
```
fitdensity(x; step, norm)
```

Obtains the density function given data sampled.

Use `step=(Float64)` to control the bin step and `norm=(0 or 1)` to determine
if the count of data points or the probability within ± step/2.

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
function fitdensity(v;nbins=nothing,
                      step=nothing, steptype="absolute",
                      vmin=nothing, vmax=nothing,
                      norm :: Int64 = 1)

  ndata = length(v)

  if vmin == nothing
    vmin = minimum(v)
  end
  if vmax == nothing
    vmax = maximum(v)
  end
  if nbins == nothing
    nbins = 100
  end
  if step == nothing
    step = (vmax - vmin)/nbins
  else
    # By default, the step size is absolute
    if steptype == "relative"
      step = step*(vmax-vmin)/nbins
    elseif steptype != "absolute"
      error(" steptype must be \"relative\" or \"absolute\"")
    end
  end

  x = Vector{Float64}(undef,nbins)
  df = Vector{Float64}(undef,nbins)

  binstep = (vmax - vmin)/nbins
  for i in 1:nbins
    x[i] = vmin + (i-1)*binstep + binstep/2
    nv = 0
    for j in 1:ndata
      if ( v[j] > x[i] - step/2 ) && ( v[j] <= x[i] + step/2 )
        nv = nv + 1
      end
    end
    binsize = min(vmax,x[i]+step/2) - max(vmin,x[i]-step/2)
    if norm == 0
      df[i] = nv
    elseif norm == 1
      df[i] = nv/(binsize*ndata)
    end
  end

  return Density(x,df,step,norm)
end
export fitdensity

function Base.show( io :: IO, d :: Density )
  if d.norm == 0
    s = "number of"
  end
  if d.norm == 1
    s = "probability of finding"
  end
  println(" 
 ------------------- Density -------------

  d contains the $s data points within x ± $(round(d.step/2,sigdigits=3))

 ----------------------------------------- ")
end
