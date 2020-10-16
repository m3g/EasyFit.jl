# 
# Linear fit
#
struct Linear
  n :: Int64
  a :: Float64
  b :: Float64
  R :: Float64
  x :: Vector{Float64}
  y :: Vector{Float64}
  resid :: Float64
  ypred :: Vector{Float64}
end

"""
`fitlinear(x,y)`

Obtains the linear fit: ``y = a*x + b``

Optional lower and upper bounds for a and b can be provided using, for example:

```
fitlinear(x,y, lower(b=0.), upper(a=5.) )
fitlinear(x,y, lower(b=0.) )
fitlinear(x,y, upper(b=0.) )
fitlinear(x,y, upper(a=1.,b=0.) )
```

# Examples
```jldoctest
julia> x = sort(rand(10)) ; y = sort(rand(10));

julia> fit = fitlinear(x,y)

------------------- Linear Fit ------------- 

Equation: y = ax + b 

With: a = 1.0448783208110997
      b = 0.18817627115683894

Pearson correlation coefficient, R = 0.8818586822210751
Average absolute residue = 0.14274752107157443

Predicted Y: ypred = [0.1987357699444139, 0.32264343301109627...

-------------------------------------------- 

```
"""
function fitlinear(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}; 
                   l :: lower = lower(), u :: upper = upper(), 
                   options :: Options = Options() )
  # Check data
  X, Y = checkdata(X,Y,options)
  # Set bounds
  vars = [ VarType(:a,Number,1), 
           VarType(:b,Number,1) ]
  lower, upper = setbounds(vars,l,u)
  # Set model
  f = (x,p) -> p[1]*x + p[2]
  ∇f = (x,p) -> [ x, 1. ]
  func = (p) -> sq_residue(p, f, X, Y)
  grad! = (p,g) -> sq_residue_grad!(p,f,∇f,X,Y,g)
  # Initial point
  p0 = Vector{Float64}(undef,2)
  initP!(p0,options,lower,upper)
  # Fit
  fit = SPGBox.spgbox!(p0, func, grad!, l=lower, u=upper)
  # Analyze results and return
  ypred = evalpred(f,X,fit)
  R = Statistics.cor(Y,ypred)
  x, y = finexy(f,X,length(X),fit)
  return Linear(length(X),fit.x...,R,x,y,fit.f,ypred)

end
@FitMethods(fitlinear)

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
  println(" Average absolute residue = ", fit.resid / fit.n)
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println("")
  println(" -------------------------------------------- ")
end

export fitlinear
