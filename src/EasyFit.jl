module EasyFit
  
  using Statistics
  using LsqFit
  using Interpolations

  Numbers = Union{Int64,Float64}
  Vectors = Union{AbstractArray{Int64},AbstractArray{Float64}}

  #
  # Random initial point
  #

  function initP(n :: Int64)
    p0 = Vector{Float64}(undef,n)
    @. p0 = -1. + 2*rand()
    return p0
  end

  # 
  # Linear fit
  #

  struct Linear
    a :: Float64
    b :: Float64
    R :: Float64
    y :: Vector{Float64}
    residue :: Vector{Float64}
  end

  function fitlinear(X :: Vectors, Y :: Vectors)
    @. model(x,p) = p[1]*x + p[2]
    p0 = initP(2)
    fit = curve_fit(model, X, Y, p0)
    y = similar(Y)
    y = model(X,fit.param)
    R = Statistics.cor(Y,y)^2
    return Linear(fit.param[1],fit.param[2],
                  R,y,fit.resid)
  end

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
    println("")
    println(" Predicted y = [",fit.y[1],", ",fit.y[2],"...")
    println(" Residues = [", fit.residue[1],", ",fit.residue[2],"...")
    println("")
    println(" -------------------------------------------- ")
  end

  export fitlinear

  #
  # Quadratic fit
  # 

  struct Quadratic
    a :: Float64
    b :: Float64
    c :: Float64
    R :: Float64
    y :: Vector{Float64}
    residue :: Vector{Float64}
  end

  function fitquadratic(X :: Vectors, Y :: Vectors)
    @. model(x,p) = p[1]*x^2 + p[2]*x + p[3]
    p0 = initP(3)
    fit = curve_fit(model, X, Y, p0)
    y = similar(Y)
    y = model(X,fit.param)
    R = Statistics.cor(Y,y)^2
    return Quadratic(fit.param[1],fit.param[2],fit.param[3],
                     R,y,fit.resid)
  end
  const fitquad = fitquadratic

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
    println(" Square Pearson correlation, R = ", fit.R)
    println("")
    println(" Predicted y = [",fit.y[1],", ",fit.y[2],"...")
    println(" Residues = [", fit.residue[1],", ",fit.residue[2],"...")
    println("")
    println(" ----------------------------------------------- ")
  end

  export fitquad, fitquadratic

  #
  # Cubic fit
  #

  struct Cubic
    a :: Float64
    b :: Float64
    c :: Float64
    d :: Float64
    R :: Float64
    y :: Vector{Float64}
    residue :: Vector{Float64}
  end

  function fitcubic(X :: Vectors, Y :: Vectors)
    @. model(x,p) = p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4]
    p0 = initP(4)
    fit = curve_fit(model, X, Y, p0)
    y = similar(Y)
    y = model(X,fit.param)
    R = Statistics.cor(Y,y)^2
    return Cubic(fit.param[1],fit.param[2],fit.param[3],fit.param[4],
                 R,y,fit.resid)
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
    println(" Square Pearson correlation, R = ", fit.R)
    println("")
    println(" Predicted y = [",fit.y[1],", ",fit.y[2],"...")
    println(" Residues = [", fit.residue[1],", ",fit.residue[2],"...")
    println("")
    println(" ----------------------------------------------- ")
  end

  export fitcubic

  #
  # Exponential fit
  #

  struct SingleExponential
    A :: Float64
    b :: Float64
    R :: Float64
    y :: Vector{Float64}
    residue :: Vector{Float64}
  end

  struct MultipleExponential
    n :: Int64
    A :: Vector{Float64}
    b :: Vector{Float64}
    R :: Float64
    y :: Vector{Float64}
    residue :: Vector{Float64}
  end

  function sum_of_exps(x :: Numbers, p :: Vector{Float64})
    n = round(Int64,length(p)/2)
    f = 0.
    for i in 1:n
      f = f + p[i]*exp(p[n+i]*x)
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

  function fitexponential(X :: Vectors, Y :: Vectors; n :: Int = 1)
    model(x,p) = exp_model(x,p)
    p0 = initP(2*n)
    fit = curve_fit(model, X, Y, p0)
    y = similar(Y)
    y = model(X,fit.param)
    R = Statistics.cor(Y,y)^2
    if n == 1
      return SingleExponential(fit.param[1],fit.param[2],
                               R,y,fit.resid)
    else
      return MultipleExponential(n,fit.param[1:n],fit.param[n+1:2*n],
                                 R,y,fit.resid)
    end
  end
  const fitexp = fitexponential

  function Base.show( io :: IO, fit :: SingleExponential )
    println("")
    println(" ------------ Single Exponential fit ----------- ")
    println("")
    println(" Equation: y = A exp(b*x) ")
    println("")
    println(" With: A = ", fit.A)
    println("       b = ", fit.b)
    println("")
    println(" Square Pearson correlation, R = ", fit.R)
    println("")
    println(" Predicted y = [",fit.y[1],", ",fit.y[2],"...")
    println(" Residues = [", fit.residue[1],", ",fit.residue[2],"...")
    println("")
    println(" ----------------------------------------------- ")
  end

  function Base.show( io :: IO, fit :: MultipleExponential )
    println("")
    println(" -------- Multiple-exponential fit ------------- ")
    println("")
    println(" Equation: y = sum(A[i] exp(b[i]) for i in 1:$(fit.n) ] ")
    println("")
    println(" With: A = ", fit.A)
    println("       b = ", fit.b)
    println("")
    println(" Square Pearson correlation, R = ", fit.R)
    println("")
    println(" Predicted y = [",fit.y[1],", ",fit.y[2],"...")
    println(" Residues = [", fit.residue[1],", ",fit.residue[2],"...")
    println("")
    println(" ----------------------------------------------- ")
  end

  export fitexp, fitexponential

  #
  # Spline
  #
 
  struct Spline
    xs :: Vector{Float64}
    ys :: Vector{Float64}
  end

  function fitspline(x :: Vectors, y :: Vectors; fine :: Int = 100)
    t = 1:length(x)
    A = hcat(x,y)
    itp = Interpolations.scale(interpolate(A, 
              (BSpline(Interpolations.Cubic(Natural(OnGrid()))), NoInterp())), t, 1:2)
    tfine = 1:length(x)/fine:length(x)
    xs, ys = [itp(t,1) for t in tfine], [itp(t,2) for t in tfine]
    return Spline(xs,ys)
  end

  function Base.show( io :: IO, fit :: Spline )
    println("")
    println(" -------- Spline fit --------------------------- ")
    println("")
    println(" x spline: xs = [",fit.xs[1],", ",fit.xs[2],"...")
    println(" y spline: ys = [",fit.ys[1],", ",fit.ys[2],"...")
    println("")
    println(" ----------------------------------------------- ")
  end

  export fitspline

end # module
