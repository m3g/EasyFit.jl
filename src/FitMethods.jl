"""

  FitMethods(f :: Func)

  Generates the methods for each fit function which allow calls using only
  upper or lower bounds to be defined, in any order.  

"""
macro FitMethods(f)
  @eval begin
    T = AbstractArray{<:Real}
    # No options provided
    $f(X :: T, Y :: T, l :: lower; s :: Float64, options :: Options) = $f(X,Y,l=l,s=s,options=options)
    $f(X :: T, Y :: T, u :: upper; kargs)= $f(X,Y,u=u,kargs...)
    $f(X :: T, Y :: T, l :: lower, u :: upper; kargs) = $f(X,Y,l=l,u=u,kargs...)
    $f(X :: T, Y :: T, u :: upper, l :: lower; kargs) = $f(X,Y,l=l,u=u,kargs...)
  end
end

macro FitMethodsExponential(f)
  @eval begin
    # No options provided
    $f(X :: T, Y :: T, l :: lower; n :: Int = 1) = $f(X,Y,l=l,n=n)
    $f(X :: T, Y :: T, u :: upper; n :: Int = 1) = $f(X,Y,u=u,n=n)
    $f(X :: T, Y :: T, l :: lower, u :: upper; n :: Int = 1) = $f(X,Y,l=l,u=u,n=n)
    $f(X :: T, Y :: T, u :: upper, l :: lower; n :: Int = 1) = $f(X,Y,l=l,u=u,n=n)
    # Options provided
    $f(X :: T, Y :: T, options :: Options; n :: Int = 1) = $f(X,Y,options=options,n=n)
    $f(X :: T, Y :: T, l :: lower, options :: Options; n :: Int = 1) = $f(X,Y,l=l,options=options,n=n)
    $f(X :: T, Y :: T, u :: upper, options :: Options; n :: Int = 1) = $f(X,Y,u=u,options=options,n=n)
    $f(X :: T, Y :: T, l :: lower, u :: upper, options :: Options; n :: Int = 1) = $f(X,Y,l=l,u=u,options=options,n=n)
    $f(X :: T, Y :: T, u :: upper, l :: lower, options :: Options; n :: Int = 1) = $f(X,Y,l=l,u=u,options=options,n=n)
  end
end
