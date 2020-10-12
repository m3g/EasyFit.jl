"""

  FitMethods(f :: Func)

  Generates the methods for each fit function which allow calls using only
  upper or lower bounds to be defined, in any order.  

"""
macro FitMethods(f)
  @eval begin
    # No options provided
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, l :: lower) = $f(X,Y,l=l)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, u :: upper) = $f(X,Y,u=u)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, l :: lower, u :: upper) = $f(X,Y,l=l,u=u)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, u :: upper, l :: lower) = $f(X,Y,l=l,u=u)
    # Options provided
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, options :: Options) = $f(X,Y,options=options)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, l :: lower, options :: Options) = $f(X,Y,l=l,options=options)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, u :: upper, options :: Options) = $f(X,Y,u=u,options=options)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, l :: lower, u :: upper, options :: Options) = $f(X,Y,l=l,u=u,options=options)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, u :: upper, l :: lower, options :: Options) = $f(X,Y,l=l,u=u,options=options)
  end
end

macro FitMethodsExponential(f)
  @eval begin
    # No options provided
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, l :: lower; n :: Int = 1) = $f(X,Y,l=l,n=n)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, u :: upper; n :: Int = 1) = $f(X,Y,u=u,n=n)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, l :: lower, u :: upper; n :: Int = 1) = $f(X,Y,l=l,u=u,n=n)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, u :: upper, l :: lower; n :: Int = 1) = $f(X,Y,l=l,u=u,n=n)
    # Options provided
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, options :: Options; n :: Int = 1) = $f(X,Y,options=options,n=n)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, l :: lower, options :: Options; n :: Int = 1) = $f(X,Y,l=l,options=options,n=n)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, u :: upper, options :: Options; n :: Int = 1) = $f(X,Y,u=u,options=options,n=n)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, l :: lower, u :: upper, options :: Options; n :: Int = 1) = $f(X,Y,l=l,u=u,options=options,n=n)
    $f(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, u :: upper, l :: lower, options :: Options; n :: Int = 1) = $f(X,Y,l=l,u=u,options=options,n=n)
  end
end
