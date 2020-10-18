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

"""
`fitquad(x,y)` or `fitquadratic(x,y)`

Obtains the quadratic fit: ``y = a*x^2 + b*x + c``

Optional lower and upper bounds for a and b can be provided using, for example:

```
fitquad(x,y, lower(b=0.), upper(a=5.,b=7.) )
fitquad(x,y, lower(b=0.) )
fitquad(x,y, upper(b=0.) )
fitquad(x,y, upper(a=1.,b=0.) )
```

# Examples
```jldoctest
julia> x = sort(rand(10)); y = rand(10).^2 .+ rand();

julia> fit = fitquad(x,y,lower(c=1))

 ------------------- Quadratic Fit ------------- 

 Equation: y = ax^2 + bx + c 

 With: a = 2.8685688648859653
       b = -1.8803961490451921
       c = 1.3269720880467206

 Pearson correlation coefficient, R = 0.5028759798430105
 Average absolute residue = 0.2137732533148255

 Predicted Y: ypred = [1.18520463048216, 1.1710010208923622...
 residues = [-0.27097835626668965, 0.2550499436259156...

 ----------------------------------------------- 

```
""" 
function fitquadratic(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}; 
                      l :: lower = lower(), u :: upper = upper(),
                      options :: Options = Options())
  # Check data
  X, Y = checkdata(X,Y,options)
  # Set bounds
  vars = [ VarType(:a,Number,1), 
           VarType(:b,Number,1),
           VarType(:c,Number,1) ]
  lower, upper = setbounds(vars,l,u)   
  # Set model
  @. model(x,p) = p[1]*x^2 + p[2]*x + p[3]
  # Fit
  fit = find_best_fit(model, X, Y, length(vars), options, lower, upper)
  # Analyze results and return
  R = pearson(X,Y,model,fit)
  x, y, ypred = finexy(X,options.fine,model,fit) 
  return Quadratic(fit.param...,R,x,y,ypred,fit.resid)
end
@FitMethods(fitquadratic)
fitquad = fitquadratic

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
  println(" Average absolute residue = ", sum(abs.(fit.residues))/length(fit.residues))
  println("")
  println(" Predicted Y: ypred = [",fit.ypred[1],", ",fit.ypred[2],"...")
  println(" residues = [", fit.residues[1],", ",fit.residues[2],"...")
  println("")
  println(" ----------------------------------------------- ")
end

export fitquad, fitquadratic
