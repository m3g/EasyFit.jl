#
# Exponential fit
#

struct SingleExponential
  ndata :: Int64
  n :: Int64
  a :: Float64
  b :: Float64
  c :: Float64
  R :: Float64
  x :: Vector{Float64}
  y :: Vector{Float64}
  resid :: Float64
  ypred :: Vector{Float64}
end

struct MultipleExponential
  ndata :: Int64
  n :: Int64
  a :: Vector{Float64}
  b :: Vector{Float64}
  c :: Float64
  R :: Float64
  x :: Vector{Float64}
  y :: Vector{Float64}
  resid :: Float64
  ypred :: Vector{Float64}
end

# Function which is a sum of exponentials
function sum_of_exps(n :: Int64, x :: Real, p :: Vector{Float64})
  np = length(p)
  f = p[np] 
  for i in 1:n
    f = f + p[i]*exp(p[n+i]*x)
  end
  return f
end

# Gradient of the sum of exponentials relative to the parameters
function grad_of_sum_of_exps(n :: Int64, x :: Real, p :: Vector{Float64})
  np = length(p)
  g = Vector{Float64}(undef,np)
  for i in 1:n
    g[i] = exp(p[n+i]*x)
    g[n+i] = x*p[i]*exp(p[n+i]*x)  
  end
  g[np] = 1.
  return g
end

"""
`fitexponential(x,y;n :: Int = 1)`  or  `fitexp(x,y;n :: Int = 1)`

Obtains single or multiexponential fits: ``y = a*exp(-x/b) + c`` or  ``y = sum(a[i]*exp(-x/b[i]) for i in 1:N) + c``

Lower and upper bounds can be optionall set:

For single exponentials all parameters are scalars:

```
fitexp(x,y,lower(c=10.),upper(b=-100))
```

For multiple exponentials, `a` and `b` bounds must be vectors of dimension `N`.

```
fitexp(x,y,n=2,lower(a=[0.,0.],c=10.),upper(b=[-100.,-5.]))
```

# Examples
```jldoctest
julia> x = sort(rand(10)); y = rand()*exp.(sort(rand(10)));

julia> fit = fitexp(x,y,lower(c=0.),n=2)

 -------- Multiple-exponential fit ------------- 

 Equation: y = sum(a[i] exp(-x/b[i]) for i in 1:2) + c 

 With: a = [-20.496181414206525, 1.2136669267088652e-7]
       b = [-45.88612863062009, 0.06979841806133893]
       c = 20.922548264186936

 Pearson correlation coefficient, R = 0.9493136213131825

 Predicted Y: ypred = [0.5925151422559018, 0.627609595906879...

 ----------------------------------------------- 

```
"""  
function fitexponential(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}; 
                        n :: Int = 1,
                        l :: lower = lower(), u :: upper = upper(),
                        options :: Options = Options())

  # Check data
  X, Y = checkdata(X,Y,options)
  # Set bounds
  if n == 1
    vars = [ VarType(:a,Number,1), 
             VarType(:b,Number,1), 
             VarType(:c,Number,1) ]
    lower, upper = setbounds(vars,l,u)     
  else
    vars = [ VarType(:a,Vector,n), 
             VarType(:b,Vector,n), 
             VarType(:c,Number,1) ]
    lower, upper = setbounds(vars,l,u)     
  end
#  for i in n+1:2*n
#    if lower[i] == -Inf
#      lower[i] = -30.
#    end
#    if upper[i] == +Inf
#      upper[i] = 30
#    end
#  end
  # Set model
#  f = (x,p) -> sum_of_exps(n,x,p)
#  ∇f = (x,p) -> grad_of_sum_of_exps(n,x,p)

  p = rand(3)
#  p = [ 1. , rand(), 0. ]

  f = (x,p) -> p[1]*exp(p[2]*x) + p[3]
  f2 = (x,p) -> log(p[1]) + p[2]*x

  ∇f1 = (x,p) -> exp(p[2]*x)
  ∇f2 = (x,p) -> x
  ∇f3 = (x,p) -> 1.

println(0," ",p)
  for i in 1:5
    newY = similar(Y)
    for i in 1:length(Y)
      newY[i] = log(Y[i] - p[3])
    end
    v = (x) -> [ p[1], x[1], p[3] ]
    fit = SPGBox.spgbox!([p[2]], 
                         x -> sq_residue(v(x),f2,X,newY), 
                         (x,g) -> sq_residue_grad!(v(x),f2,∇f2,X,newY,g),
                         l=[lower[2]], u=[upper[2]])
    p[2] = fit.x[1]
println(2," ",p," ",fit.f)

    v = (x) -> [ x[1], p[2], p[3] ]
    fit = SPGBox.spgbox!([p[1]], 
                         x -> sq_residue(v(x),f,X,Y), 
                         (x,g) -> sq_residue_grad!(v(x),f,∇f1,X,Y,g), 
                         l=[lower[1]], u=[upper[1]])
    p[1] = fit.x[1]
println(1," ",p," ",fit.f)


    v = (x) -> [ p[1], p[2], x[1] ]
    fit = SPGBox.spgbox!([p[3]], 
                         x -> sq_residue(v(x),f,X,Y), 
                         (x,g) -> sq_residue_grad!(v(x),f,∇f3,X,Y,g),
                         l=[lower[3]], u=[upper[3]])
    p[3] = fit.x[1]
println(3," ",p," ",fit.f)
  end




error("FIM")


  # Fit
  fit = find_best_fit(f, ∇f, X, Y, length(vars), options, lower, upper)
  # Analyze results and return
  ypred = evalpred(f,X,fit)
  R = Statistics.cor(Y,ypred)
  x, y = finexy(f,X,length(X),fit)   
  a = fit.x[1:n]
  b = fit.x[n+1:2*n]
  c = fit.x[2*n+1]
  if n == 1
    return SingleExponential(length(X),1,a[1],b[1],c,R,x,y,fit.f,ypred)
  else
    return MultipleExponential(length(X),n,a,b,c,R,x,y,fit.f,ypred)
  end
end
@FitMethodsExponential(fitexponential)
fitexp = fitexponential

function Base.show( io :: IO, fit :: SingleExponential )
  println("")
  println(" ------------ Single Exponential fit ----------- ")
  println("")
  println(" Equation: y = a × exp(-x/b) + c")
  println("")
  println(" With: a = ", fit.a)
  println("       b = ", fit.b)
  println("       c = ", fit.c)
  println("")
  println(" Pearson correlation coefficient, R = ", fit.R)
  println(" Average absolute square residue = ", fit.resid / fit.ndata)
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println("")
  println(" ----------------------------------------------- ")
end

function Base.show( io :: IO, fit :: MultipleExponential )
  println("")
  println(" -------- Multiple-exponential fit ------------- ")
  println("")
  println(" Equation: y = sum(a[i] × exp(-x/b[i]) for i in 1:$(fit.n)) + c ")
  println("")
  println(" With: a = ", fit.a)
  println("       b = ", fit.b)
  println("       c = ", fit.c)
  println("")
  println(" Pearson correlation coefficient, R = ", fit.R)
  println(" Average absolute square residue = ", fit.resid / fit.ndata)
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println("")
  println(" ----------------------------------------------- ")
end

export fitexp, fitexponential
