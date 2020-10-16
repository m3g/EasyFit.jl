#
# Cubic fit
#

struct Cubic
  n :: Int64
  a :: Float64
  b :: Float64
  c :: Float64
  d :: Float64
  R :: Float64
  x :: Vector{Float64}
  y :: Vector{Float64}
  resid :: Float64
  ypred :: Vector{Float64}
end

"""
`fitcubic(x,y)`

Obtains the cubic polynomial fit: ``y = ax^3 + bx^2 + cx + d``

Optional lower and upper bounds for a and b can be provided using, for example:

```
fitcubic(x,y, lower(b=0.), upper(a=5.,d=10.) )
fitcubic(x,y, lower(b=0.,c=0.) )
fitcubic(x,y, upper(b=0.) )
fitcubic(x,y, upper(a=1.,b=0.) )
```

# Examples
```jldoctest
julia> x = sort(rand(10)) ; y = sort(rand(10)).^3;

julia> fit = fitcubic(x,y,lower(d=5.),upper(d=6.))

 ------------------- Cubic Fit ----------------- 

 Equation: y = ax^3 + bx^2 + cx + d 

 With: a = 4.797115186839406
       b = -5.268081309993556
       c = 2.0751456266896042
       d = -0.042233689646186394

 Pearson correlation coefficient, R = 0.9900635473581005
 Average absolute residue = 0.0352610880828611

 Predicted Y: ypred = [-0.03381603021235922, 0.026984191488563257...

 ----------------------------------------------- 

```
"""   
function fitcubic(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}; 
                  l :: lower = lower(), u :: upper = upper(),
                  options :: Options = Options())
  # Check data
  X, Y = checkdata(X,Y,options)
  # Set bounds
  vars = [ VarType(:a,Number,1),
           VarType(:b,Number,1),
           VarType(:c,Number,1),
           VarType(:d,Number,1) ]
  lower, upper = setbounds(vars,l,u)
  # Set model
  f = (x,p) -> p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4]
  ∇f = (x,p) -> [ x^3, x^2, x, 1. ]
  # Fit
  fit = find_best_fit(f, ∇f, X, Y, length(vars), options, lower, upper)
  # Analyze results and return
  ypred = evalpred(f,X,fit)
  R = Statistics.cor(Y,ypred)
  x, y = finexy(f,X,length(X),fit)  
  return Cubic(length(X),fit.x...,R,x,y,fit.f,ypred)    
end
@FitMethods(fitcubic)

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
  println(" Average absolute residue = ", fit.resid/fit.n)
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println("")
  println(" ----------------------------------------------- ")
end

export fitcubic
