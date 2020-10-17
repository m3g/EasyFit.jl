#
# Computes the moving average of the data
#

struct MovingAverage
  n :: Int64
  x :: Vector{Float64}
  R :: Float64
  residues :: Vector{Float64}
end

"""

`movingaverage(x,n)` or `movavg(x,n)`

Computes a moving average of `x[i]` in the range `i±(n-1)/2`. If `n` is even, we do `n ←  n + 1`.

# Example

```jldoctest
julia> x = rand(10);

julia> movingaverage(x,3)

 ------------------- Moving Average ----------

 Number of points averaged: 3 (± 1 points)

 Pearson correlation coefficient, R = 0.3532754137104625

 Averaged X: x = [0.5807828672543551, 0.40496733381946143...
 residues = [-0.22791917753944557, 0.4037347109743393...

 --------------------------------------------
```

"""
function movingaverage( X :: AbstractArray{<:Real}, n :: Int )
  if n == 0
    error(" Please set n with movavg(x,n=10), for example. ")
  end
  if ! isodd(n)
    n = n + 1
  end
  delta = round(Int64,(n-1)/2)
  nX = length(X)
  x = similar(X) 
  for i in 1:nX
    jmin = max(i-delta,1)
    jmax = min(i+delta,nX)
    x[i] = 0.
    for j in jmin:jmax
      x[i] = x[i] + X[j]
    end
    x[i] = x[i] / (jmax-jmin+1)
  end
  residues = similar(X)
  for i in 1:nX
    residues[i] = X[i] - x[i]
  end
  R = Statistics.cor(x,X)
  return MovingAverage(n,x,R,residues)
end
movingaverage(X::AbstractArray{<:Real};n=Int=0) = movingaverage(X,n)
movavg = movingaverage

function Base.show( io :: IO, fit :: MovingAverage )
  println("")
  println(" ------------------- Moving Average ---------- ")
  println("")
  println(" Number of points averaged: $(fit.n) (± $(round(Int64,(fit.n-1)/2)) points)")
  println("")
  println(" Pearson correlation coefficient, R = ", fit.R)
  println("")
  println(" Averaged X: x = [",fit.x[1],", ",fit.x[2],"...")
  println(" residues = [", fit.residues[1],", ",fit.residues[2],"...")
  println("")
  println(" -------------------------------------------- ")
end

export movingaverage, movavg

