#
# Cubic fit
#

struct Cubic
  a :: Float64
  b :: Float64
  c :: Float64
  d :: Float64
  R :: Float64
  x :: Vector{Float64}
  y :: Vector{Float64}
  ypred :: Vector{Float64}
  residues :: Vector{Float64}
end

function fitcubic(X :: Vectors, Y :: Vectors, options :: Options)
  @. model(x,p) = p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4]
  fit = find_best_fit(model, X, Y, 4, options)
  R = pearson(X,Y,model,fit)
  x, y, ypred = finexy(X,options.fine,model,fit) 
  return Cubic(fit.param...,R,x,y,ypred,fit.resid)
end
fitcubic(X :: Vectors, Y :: Vectors) = fitcubic(X,Y,Options())

function Base.show( io :: IO, fit :: Cubic )
  println("")
  println(" ------------------- Cubic Fit ----------------- ")
  println("")
  println(" Equation: y = ax^3 + bx^2 + cx + d ")
  println("")
  println(" With: a = ", fit.a)
  println("       b = ", fit.b)
  println("       c = ", fit.c)
  println("       d = ", fit.d)
  println("")
  println(" Pearson correlation coefficient, R = ", fit.R)
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println(" residues = [", fit.residues[1],", ",fit.residues[2],"...")
  println("")
  println(" ----------------------------------------------- ")
end

export fitcubic
