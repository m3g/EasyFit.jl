#
# Computes the moving average of the data
#

struct MovingAverage
  n :: Int64
  x :: Vector{Float64}
  y :: Vector{Float64}
  R :: Float64
  residues :: Vector{Float64}
end

function movingaverage( X :: Vectors, Y :: Vectors, n :: Int )
  if ! isodd(n)
    n = n + 1
  end
  delta = round(Int64,(n-1)/2)
  nY = length(Y)
  y = similar(Y) 
  for i in 1:nY
    jmin = max(i-delta,1)
    jmax = min(i+delta,nY)
    y[i] = 0.
    for j in jmin:jmax
      y[i] = y[i] + Y[j]
    end
    y[i] = y[i] / (jmax-jmin+1)
  end
  residues = similar(Y)
  for i in 1:nY
    residues[i] = Y[i] - y[i]
  end
  R = Statistics.cor(y,Y)
  return MovingAverage(n,X,y,R,residues)
end
function movingaverage(X :: Vectors, Y :: Vectors; n :: Int = 0)
  if n == 0
    error(" Please set n with movavg(x,y,n=10), for example. ")
  end
  return movingaverage(X,Y,n)
end
movavg(X :: Vectors, Y :: Vectors, n :: Int) = movingaverage(X,Y,n)
movavg(X :: Vectors, Y :: Vectors; n :: Int = 0) = movingaverage(X,Y,n) 

function Base.show( io :: IO, fit :: MovingAverage )
  println("")
  println(" ------------------- Moving Average ---------- ")
  println("")
  println(" Number of points averaged: $(fit.n) (± $(round(Int64,(fit.n-1)/2)) points)")
  println("")
  println(" Pearson correlation coefficient, R = ", fit.R)
  println("")
  println(" Averaged Y: y = [",fit.y[1],", ",fit.y[2],"...")
  println(" residues = [", fit.residues[1],", ",fit.residues[2],"...")
  println("")
  println(" -------------------------------------------- ")
end

export movingaverage, movavg

