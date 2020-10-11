#
# Exponential fit
#

struct SingleExponential
  a :: Float64
  b :: Float64
  c :: Float64
  R :: Float64
  x :: Vector{Float64}
  y :: Vector{Float64}
  ypred :: Vector{Float64}
  residues :: Vector{Float64}
end

struct MultipleExponential
  n :: Int64
  a :: Vector{Float64}
  b :: Vector{Float64}
  c :: Float64
  R :: Float64
  x :: Vector{Float64}
  y :: Vector{Float64}
  ypred :: Vector{Float64}
  residues :: Vector{Float64}
end

function sum_of_exps(x :: Real, p :: Vector{Float64})
  n = round(Int64,(length(p)-1)/2)
  f = p[length(p)] 
  for i in 1:n
    f = f + p[i]*exp(x/p[n+i])
  end
  return f
end

function exp_model(X :: AbstractArray{<:Real}, p :: Vector{Float64}) 
  f = Vector{Float64}(undef,length(X))  
  for i in 1:length(X)
    f[i] = sum_of_exps(X[i],p)
  end
  return f
end

"""
`fitexponential(x,y;n :: Int = 1)`  or  `fitexp(x,y;n :: Int = 1)`

Obtains single or multiexponential fits: ``y = a*exp(x/b) + c`` or  ``y = sum(a[i]*exp(x/b[i]) for i in 1:N) + c``

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

julia> fitexp(x,y,lower(c=0.),n=2)

 -------- Multiple-exponential fit ------------- 

 Equation: y = sum(a[i] exp(x/b[i]) for i in 1:2) + c 

 With: a = [-20.496181414206525, 1.2136669267088652e-7]
       b = [-45.88612863062009, 0.06979841806133893]
       c = 20.922548264186936

 Pearson correlation coefficient, R = 0.9493136213131825

 Predicted Y: ypred = [0.5925151422559018, 0.627609595906879...
 residues = [0.04784868767934802, 0.02061037608043914...

 ----------------------------------------------- 

```
"""  
function fitexponential(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}; 
                        n :: Int = 1,
                        l :: lower = lower(), u :: upper = upper(),
                        options :: Options = Options())

  # Check data
  X, Y = checkdata(X,Y,options)
  # Set model
  model(x,p) = exp_model(x,p)
  # Number of exponentials
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
    a = fit.param[1:n][ind]
    b = fit.param[n+1:2*n][ind]
    return MultipleExponential(n,a,b,fit.param[2*n+1],
                               R,x,y,ypred,fit.resid)
  end
end
@FitMethodsExponential(fitexponential)
fitexp = fitexponential

function Base.show( io :: IO, fit :: SingleExponential )
  println("")
  println(" ------------ Single Exponential fit ----------- ")
  println("")
  println(" Equation: y = a exp(x/b) + c")
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

function Base.show( io :: IO, fit :: MultipleExponential )
  println("")
  println(" -------- Multiple-exponential fit ------------- ")
  println("")
  println(" Equation: y = sum(a[i] exp(x/b[i]) for i in 1:$(fit.n)) + c ")
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

export fitexp, fitexponential
