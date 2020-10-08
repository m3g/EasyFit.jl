# 
# Linear fit
#

struct Linear
  a :: Float64
  b :: Float64
  R :: Float64
  x :: Vector{Float64}
  y :: Vector{Float64}
  ypred :: Vector{Float64}
  residues :: Vector{Float64}
end

function fitlinear(X :: Vectors, Y :: Vectors, options :: Options)
  X, Y = checkdata(X,Y)
  @. model(x,p) = p[1]*x + p[2]
  p0 = initP(2,options)
  fit = curve_fit(model, X, Y, p0)
  R = pearson(X,Y,model,fit)
  x, y, ypred = finexy(X,length(X),model,fit)
  return Linear(fit.param...,R,x,y,ypred,fit.resid)
end
fitlinear(X :: Vectors, Y :: Vectors) = fitlinear(X,Y,Options())

function Base.show( io :: IO, fit :: Linear )
  println("")
  println(" ------------------- Linear Fit ------------- ")
  println("")
  println(" Equation: y = ax + b ")
  println("")
  println(" With: a = ", fit.a)
  println("       b = ", fit.b)
  println("")
  println(" Pearson correlation coefficient, R = ", fit.R)
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println(" residues = [", fit.residues[1],", ",fit.residues[2],"...")
  println("")
  println(" -------------------------------------------- ")
end

export fitlinear
