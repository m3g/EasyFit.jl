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

"""
`fitcubic(x,y)`

Obtains the cubic polynomial fit: ``y = ax^3 + bx^2 + cx + d``

Optional lower and upper bounds for a, b, and c can be provided using, for example:

```
fitcubic(x,y, l=lower(b=0.), u=upper(a=5.))
```

and `d` can be set to constant with, for example:

```
fitcubic(x,y,d=5.)
```

# Examples
```jldoctest
julia>  x = sort(rand(10)); y = x.^3 .+ rand(10);

julia> fit = fitcubic(x,y)

 ------------------- Cubic Fit ----------------- 

 Equation: y = ax^3 + bx^2 + cx + d 

 With: a = 12.637633791600711
       b = -19.648194970330454
       c = 10.018385827387148
       d = -0.8740912356800155

 Pearson correlation coefficient, R = 0.7831345513024988
 Average square residue = 0.0781543071776559

 Predicted Y: ypred = [0.24999805379642903, 0.3001612840610868...
 residues = [0.2238223147726266, 0.12656861200050698...

 ----------------------------------------------- 

```
"""   
function fitcubic(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}; 
                  l :: lower = lower(), u :: upper = upper(), d = nothing,
                  options :: Options = Options())
  # Check data
  X, Y = checkdata(X,Y,options)
  # Set bounds
  vars = [ VarType(:a,Number,1),
           VarType(:b,Number,1),
           VarType(:c,Number,1),
           VarType(:d,Nothing,1) ]
  lower, upper = setbounds(vars,l,u)
  if d == nothing
    # Set model
    @. model(x,p) = p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4]
    # Fit
    fit = find_best_fit(model, X, Y, length(vars), options, lower, upper)
    # Analyze results and return
    R = pearson(X,Y,model,fit)
    x, y, ypred = finexy(X,options.fine,model,fit) 
    return Cubic(fit.param...,R,x,y,ypred,fit.resid)
  else
    lower = lower[1:length(vars)-1]
    upper = upper[1:length(vars)-1]
    # Set model
    @. model_const(x,p) = p[1]*x^3 + p[2]*x^2 + p[3]*x + d
    # Fit
    fit = find_best_fit(model_const, X, Y, length(vars)-1, options, lower, upper)
    # Analyze results and return
    R = pearson(X,Y,model_const,fit)
    x, y, ypred = finexy(X,options.fine,model,fit) 
    return Cubic(fit.param...,d,R,x,y,ypred,fit.resid)
  end
end

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
  println(" Average square residue = ", mean(fit.residues.^2))
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println(" residues = [", fit.residues[1],", ",fit.residues[2],"...")
  println("")
  println(" ----------------------------------------------- ")
end

export fitcubic
