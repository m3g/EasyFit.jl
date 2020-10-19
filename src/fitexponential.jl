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
  n = div(length(p)-1,2)
  f = p[length(p)] 
  for i in 1:n
    f = f + p[i]*exp(-x/p[n+i])
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

function sum_of_exps_const(x :: Real, p :: Vector{Float64},c)
  n = div(length(p),2)
  f = c
  for i in 1:n
    f = f + p[i]*exp(-x/p[n+i])
  end
  return f
end

function exp_model_const(X :: AbstractArray{<:Real}, p :: Vector{Float64},c) 
  f = Vector{Float64}(undef,length(X))  
  for i in 1:length(X)
    f[i] = sum_of_exps_const(X[i],p,c)
  end
  return f
end

"""
`fitexponential(x,y;n :: Int = 1)`  or  `fitexp(x,y;n :: Int = 1)`

Obtains single or multiexponential fits: ``y = a*exp(-x/b) + c`` or  ``y = sum(a[i]*exp(-x/b[i]) for i in 1:N) + c``

Lower and upper bounds can be optionall set, and the intercept `c` can set to be constant:

For single exponentials all parameters are scalars:

```
fitexp(x,y,l=lower(b=0.),u=upper(b=10.),c=5.)
```

For multiple exponentials, `a` and `b` bounds must be vectors of dimension `N`.

```
fitexp(x,y,n=2,l=lower(a=[0.,0.]),u=upper(b=[-100.,-5.]))
```

# Examples
```jldoctest
julia> x = sort(rand(10)); y = rand()*exp.(sort(rand(10)));

julia> fit = fitexp(x,y,l=lower(b=[0.,0.]),n=2)

 -------- Multiple-exponential fit -------------

 Equation: y = sum(a[i] exp(-x/b[i]) for i in 1:2) + c

 With: a = [6.60693727987886e-13, 0.6249999999993409]
       b = [0.02688289803014393, 0.5000000000002596]
       c = 0.37499999999999856

 Pearson correlation coefficient, R = 1.0
 Average square residue = 1.1639900380979497e-29

 Predicted Y: ypred = [1.0000000000000002, 0.4595845520228801...
 residues = [2.220446049250313e-16, -2.831068712794149e-15...

 -----------------------------------------------

```
"""  
function fitexponential(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}; 
                        n :: Int = 1, c = nothing,
                        l :: lower = lower(), u :: upper = upper(),
                        options :: Options = Options())

  # Check data
  X, Y = checkdata(X,Y,options)
  if c == nothing
    # Set model
    model(x,p) = exp_model(x,p)
    # Number of exponentials
    # Set bounds
    if n == 1
      vars = [ VarType(:a,Number,1), 
               VarType(:b,Number,1), 
               VarType(:c,Nothing,1) ]
      lower, upper = setbounds(vars,l,u)     
    else
      vars = [ VarType(:a,Vector,n), 
               VarType(:b,Vector,n), 
               VarType(:c,Nothing,1) ]
      lower, upper = setbounds(vars,l,u)     
    end
    # Fit
    fit = find_best_fit(model,X,Y,2*n+1,options,lower,upper)
    # Analyze and return
    R = pearson(X,Y,model,fit)
    x, y, ypred = finexy(X,options.fine,model,fit) 
    if n == 1
      return SingleExponential(length(X),fit.param[1],fit.param[2],fit.param[3],
                               R,x,y,ypred,fit.resid)
    else
      ind = collect(1:n)
      sort!(ind,by=i->fit.param[n+i])
      a = fit.param[1:n][ind]
      b = fit.param[n+1:2*n][ind]
      return MultipleExponential(length(X),n,a,b,fit.param[2*n+1],
                                 R,x,y,ypred,fit.resid)
    end
  else
    # Set model
    model_const(x,p) = exp_model_const(x,p,c)
    # Number of exponentials
    # Set bounds
    if n == 1
      vars = [ VarType(:a,Number,1), 
               VarType(:b,Number,1), 
               VarType(:c,Nothing,1) ]
      lower, upper = setbounds(vars,l,u)     
    else
      vars = [ VarType(:a,Vector,n), 
               VarType(:b,Vector,n), 
               VarType(:c,Nothing,1) ]
      lower, upper = setbounds(vars,l,u)     
    end
    lower = lower[1:2*n]
    upper = upper[1:2*n]
    # Fit
    fit = find_best_fit(model_const,X,Y,2*n,options,lower,upper)
    # Analyze and return
    R = pearson(X,Y,model_const,fit)
    x, y, ypred = finexy(X,options.fine,model_const,fit) 
    if n == 1
      return SingleExponential(fit.param[1],fit.param[2],c,R,x,y,ypred,fit.resid)
    else
      ind = collect(1:n)
      sort!(ind,by=i->fit.param[n+i])
      a = fit.param[1:n][ind]
      b = fit.param[n+1:2*n][ind]
      return MultipleExponential(n,a,b,c,R,x,y,ypred,fit.resid)
    end
  end
end
fitexp = fitexponential

function Base.show( io :: IO, fit :: SingleExponential )
  println("")
  println(" ------------ Single Exponential fit ----------- ")
  println("")
  println(" Equation: y = a exp(-x/b) + c")
  println("")
  println(" With: a = ", fit.a)
  println("       b = ", fit.b)
  println("       c = ", fit.c)
  println("")
  println(" Pearson correlation coefficient, R = ", fit.R)
  println(" Average square residue = ", mean(fit.residues.^2))
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
  println(" Equation: y = sum(a[i] exp(-x/b[i]) for i in 1:$(fit.n)) + c ")
  println("")
  println(" With: a = ", fit.a)
  println("       b = ", fit.b)
  println("       c = ", fit.c)
  println("")
  println(" Pearson correlation coefficient, R = ", fit.R)
  println(" Average square residue = ", mean(fit.residues.^2))
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println(" residues = [", fit.residues[1],", ",fit.residues[2],"...")
  println("")
  println(" ----------------------------------------------- ")
end

export fitexp, fitexponential
