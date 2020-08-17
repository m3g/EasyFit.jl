# EasyFit

Easy interface for obtaining fits for 2D data.

The purpose of this package is to provide a very simple interface to obtain
some of the most common fits to 2D data. Currently, simple fitting functions
are available for linear, quadratic, cubic, exponential, and spline fits.

On the background this interface uses the [LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl)
and [Interpolations](https://github.com/JuliaMath/Interpolations.jl), which are already
quite easy to use.

Our aim is to provide a package for quick fits without having to think about the code.

## Installation

```
julia> ] add EasyFit

julia> using EasyFit

```

## Contents

- [Linear fit](#linear)

<a name="linear"/>

## Linear fit

To perform a linear fitting, do:

```julia
julia> x = sort(rand(10)); y = sort(rand(10)); # some data

julia> fit = fitlinear(x,y)

 ------------------- Linear Fit ------------- 

 Equation: y = ax + b 

 With: a = 0.605935046984355
       b = -0.027026434622431462

 Pearson correlation coefficient, R = 0.5630137802184002

 Predicted y = [-0.009488291459872424, -0.004421217036880542...
 Residues = [-0.08666886144188027, -0.12771486789744962...

 -------------------------------------------- 

```

The `fit` data structure which comes out of `fitlinear` contains the output data with
the same names as shown in the above output:

```julia
julia> fit.a
0.9245529646308137

julia> fit.b
0.08608398402393584

julia> fit.R
0.765338307304594

julia> fit.y[1:3]
3-element Array{Float64,1}:
 0.14440450322996096
 0.18487629923663051
 0.27846925434751973

julia> fit.residue[1:3]
3-element Array{Float64,1}:
 0.143300680784536
 0.16615578533670206
 0.15671338348944602


```




*Quadratic fitting*

```julia
julia> x = sort(rand(10)); y = sort(rand(10)).^2; # some data

julia> fitquad(x,y)  # or fitquadratic(x,y)

 ------------------- Quadratic Fit ------------- 

 Equation: y = ax^2 + bx + c 

 With: a = 0.935408728589832
       b = 0.07985866671623199
       c = 0.08681962205579699

 Square Pearson correlation, R = 0.9591045325089623

 Predicted y = [0.08910633345247763, 0.08943732276526263...
 Residues = [0.07660191693062993, 0.07143385689027287...

 ----------------------------------------------- 

```

*Other fitting functions*

Other fitting functions avaiable are:

```julia
julia> fitcubic(x,y)
          ...
julia> fitexp(x,y) # or fitexponential(x,y)
          ...
julia> fitexp(x,y,n=3) # or fitexponential(x,y,n=3) -- for multiple exponentials 
          ...
julia> fitspline(x,y) # smoothness can be controled by fitspline(x,y,fine=1000)
          ...

```



