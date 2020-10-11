#
# Exponential fit
#

struct SingleExponential
  A :: Float64
  b :: Float64
  C :: Float64
  R :: Float64
  x :: Vector{Float64}
  y :: Vector{Float64}
  ypred :: Vector{Float64}
  residues :: Vector{Float64}
end

struct MultipleExponential
  n :: Int64
  A :: Vector{Float64}
  b :: Vector{Float64}
  C :: Float64
  R :: Float64
  x :: Vector{Float64}
  y :: Vector{Float64}
  ypred :: Vector{Float64}
  residues :: Vector{Float64}
end

function sum_of_exps(x :: Numbers, p :: Vector{Float64})
  n = round(Int64,(length(p)-1)/2)
  f = p[length(p)] 
  for i in 1:n
    f = f + p[i]*exp(x/p[n+i])
  end
  return f
end

function exp_model(X :: Vectors, p :: Vector{Float64}) 
  f = Vector{Float64}(undef,length(X))  
  for i in 1:length(X)
    f[i] = sum_of_exps(X[i],p)
  end
  return f
end

function fitexponential(X :: Vectors, Y :: Vectors, options :: Options; n :: Int = 1)
  # Check data
  X, Y = checkdata(X,Y)
  # Check bounds
  lower = Vector{Float64}(undef,2*n+1)
  upper = Vector{Float64}(undef,2*n+1)
  if n == 1
    lower[1] = lower_bound(options.lower.A,"A",Number)
    lower[2] = lower_bound(options.lower.b,"b",Number)
    lower[3] = lower_bound(options.lower.C,"C",Number)
    upper[1] = upper_bound(options.upper.A,"A",Number)
    upper[2] = upper_bound(options.upper.b,"b",Number)
    upper[3] = upper_bound(options.upper.C,"C",Number)
  else
    @. lower[1:n] = lower_bound(options.lower.A,"A",Vector)
    @. lower[n+1:2*n] = lower_bound(options.lower.b,"b",Vector)
    lower[2*n+1] = lower_bound(options.lower.C,"C",Number)
    @. upper[1:n] = upper_bound(options.upper.A,"A",Vector)
    @. upper[n+1:2*n] = upper_bound(options.upper.b,"b",Vector)
    upper[2*n+1] = upper_bound(options.upper.C,"C",Number)
  end
  # Model
  model(x,p) = exp_model(x,p)
  # Fit
  fit = find_best_fit(model,X,Y,2*n+1,options,lower,upper)
  # Analyze and return
  R = pearson(X,Y,model,fit)
  x, y, ypred = finexy(X,options.fine,model,fit) 
  if n == 1
    return SingleExponential(fit.param[1],fit.param[2],fit.param[3],
                             R,x,y,ypred,fit.resid)
  else
    ind = collect(1:n)
    sort!(ind,by=i->fit.param[n+i])
    A = fit.param[1:n][ind]
    b = fit.param[n+1:2*n][ind]
    return MultipleExponential(n,A,b,fit.param[2*n+1],
                               R,x,y,ypred,fit.resid)
  end
end
fitexponential(X :: Vectors, Y :: Vectors; n :: Int = 1) = fitexponential(X,Y,Options(),n=n)
fitexp(X :: Vectors, Y :: Vectors; n :: Int = 1) = fitexponential(X,Y,Options(),n=n)
fitexp(X :: Vectors, Y :: Vectors, options :: Options; n :: Int = 1) = fitexponential(X,Y,options,n=n)

function Base.show( io :: IO, fit :: SingleExponential )
  println("")
  println(" ------------ Single Exponential fit ----------- ")
  println("")
  println(" Equation: y = A exp(x/b) + C")
  println("")
  println(" With: A = ", fit.A)
  println("       b = ", fit.b)
  println("       C = ", fit.C)
  println("")
  println(" Pearson correlation coefficient, R = ", fit.R)
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println(" residues = [", fit.residues[1],", ",fit.residues[2],"...")
  println("")
  println(" ----------------------------------------------- ")
end

function Base.show( io :: IO, fit :: MultipleExponential )
  println("")
  println(" -------- Multiple-exponential fit ------------- ")
  println("")
  println(" Equation: y = sum(A[i] exp(x/b[i]) for i in 1:$(fit.n)) + C ")
  println("")
  println(" With: A = ", fit.A)
  println("       b = ", fit.b)
  println("       C = ", fit.C)
  println("")
  println(" Pearson correlation coefficient, R = ", fit.R)
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println(" residues = [", fit.residues[1],", ",fit.residues[2],"...")
  println("")
  println(" ----------------------------------------------- ")
end

function fitexp() 
  println(" if n = 1: Equation: y = A exp(x/b) + C")
  println(" if n > 1: Equation: y = sum(A[i] exp(x/b[i]) for i in 1:n + C ")
end

export fitexp, fitexponential
