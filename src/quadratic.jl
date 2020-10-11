#
# Quadratic fit
# 

struct Quadratic
  a :: Float64
  b :: Float64
  c :: Float64
  R :: Float64
  x :: Vector{Float64}
  y :: Vector{Float64}
  ypred :: Vector{Float64}
  residues :: Vector{Float64}
end

function fitquadratic(X :: Vectors, Y :: Vectors, options :: Options)
  X, Y = checkdata(X,Y)
  @. model(x,p) = p[1]*x^2 + p[2]*x + p[3]
  fit = find_best_fit(model, X, Y, 3, options)
  R = pearson(X,Y,model,fit)
  x, y, ypred = finexy(X,options.fine,model,fit) 
  return Quadratic(fit.param...,R,x,y,ypred,fit.resid)
end
fitquadratic(X :: Vectors, Y :: Vectors) = fitquadratic(X,Y,Options())
fitquad(X :: Vectors, Y :: Vectors) = fitquadratic(X,Y,Options())
fitquad(X :: Vectors, Y :: Vectors, options :: Options) = fitquadratic(X,Y,options)

function Base.show( io :: IO, fit :: Quadratic )
  println("")
  println(" ------------------- Quadratic Fit ------------- ")
  println("")
  println(" Equation: y = ax^2 + bx + c ")
  println("")
  println(" With: a = ", fit.a)
  println("       b = ", fit.b)
  println("       c = ", fit.c)
  println("")
  println(" Pearson correlation coefficient, R = ", fit.R)
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println(" residues = [", fit.residues[1],", ",fit.residues[2],"...")
  println("")
  println(" ----------------------------------------------- ")
end

fitquadratic() = println(" Equation: y = ax^2 + bx + c ")
fitquad() = fitquadratic()

export fitquad, fitquadratic
