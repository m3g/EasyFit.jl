#
# Spline
#

struct Spline
  x :: Vector{Float64}
  y :: Vector{Float64}
end

function fitspline(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, options :: Options)
  X, Y = checkdata(X,Y)
  t = 1:length(X)
  A = hcat(X,Y)
  itp = Interpolations.scale(interpolate(A, 
            (BSpline(Interpolations.Cubic(Natural(OnGrid()))), NoInterp())), t, 1:2)
  tfine = 1:length(X)/options.fine:length(X)
  x, y = [itp(t,1) for t in tfine], [itp(t,2) for t in tfine]
  return Spline(x,y)
end
fitspline(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}) = fitspline(X,Y,Options())

function Base.show( io :: IO, fit :: Spline )
  println("")
  println(" -------- Spline fit --------------------------- ")
  println("")
  println(" x spline: x = [",fit.x[1],", ",fit.x[2],"...")
  println(" y spline: y = [",fit.y[1],", ",fit.y[2],"...")
  println("")
  println(" ----------------------------------------------- ")
end

export fitspline
