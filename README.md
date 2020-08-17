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

Read the `Linear Fit` section first, because all the others are similar, with
few specificities:

- [Example output](#example)
- [Linear fit](#linear)
- [Quadratic fit](#quad)
- [Cubic fit](#cubic)
- [Exponential fit](#exp)
- [Splines](#splines)

All functions except the linear fits accept an additional keyword called `fine`, which 
determines how many points the output vectors of the fit will have, such that the plots
of the fits are smooth. By default, `fine=100 , which means that the fits will have 
100 times the number of points of the original data. 

<a name="linear"/>

## Example output:

This figure was obtained using `Plots`, after obtaining a fit of each type, with

```julia
julia> scatter!(x,y) # plot original data
julia> plot!(fit.x,fit.y) # plot the resulting fit

```

<img src="https://raw.githubusercontent.com/m3g/EasyFit/master/docs/plots.png">

The complete script is available at: 
[plots.jl](#https://raw.githubusercontent.com/m3g/EasyFit/master/docs/plots.jl)


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


```

The `fit.x` and `fit.y` vectors can be used for plotting the results:

```julia
julia> using Plots

julia> scatter(x,y) # the original data

julia> plot!(fit.x,fit.y) # the fit


```

<a name="quad"/>

## Quadratic fit

Use the `fitquad` function:

```julia
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

<a name="cubic"/>

## Cubic fit

Use the `fitcubic` function:

```julia
julia> fitcubic(x,y) 

 ------------------- Cubic Fit ----------------- 

 Equation: y = ax^3 + bx^2 + cx + d 

 With: a = 1.6860182468269271
       b = -2.197790204605215
       c = 1.431666717127516
       d = -0.10389199522825227

 Square Pearson correlation, R = 0.9636627003609293

 Predicted Y: ypred = [0.024757602237563042, 0.1762724543346461...
 residues = [-0.021614675301217884, 0.0668157091306878...

 ----------------------------------------------- 


```

<a name="exp"/>

## Exponential fits

Use the `fitexp` function:

```julia
julia> fitexp(x,y) # or fitexponential

 ------------ Single Exponential fit ----------- 

 Equation: y = A exp(x/b) 

 With: A = 0.08309782657193134
       b = 0.4408664103095801

 Square Pearson correlation, R = 0.957162526367073

 Predicted Y: ypred = [0.10558554154948542, 0.16605481935145136...
 residues = [0.059213264010704494, 0.056598074147493044...

 ----------------------------------------------- 


```

or add `n=N` for a multiple-exponential fit:

```julia

julia> fit = EasyFit.fitexp(x,y,n=3)

 -------- Multiple-exponential fit ------------- 

 Equation: y = sum(A[i] exp(x/b[i]) for i in 1:3 ] 

 With: A = [2.0447736471832363e-16, 3.165225832379937, -3.2171314371600785]
       b = [0.02763465220057311, -46969.25088088338, -4.403370258345724]

 Square Pearson correlation, R = 0.9835776339692254

 Predicted Y: ypred = [0.024313571992034433, 0.1635108558614995...
 residues = [-0.022058705546746493, 0.05405411065754118...

 ----------------------------------------------- 

```

!!! Warning: exponential fits can be tricky. Run multiple times (which
    generates new initial points) if you don't like the result.


<a name="spline"/>

## Splines

Use the `fitspline` function:

```julia

julia> fit = EasyFit.fitspline(x,y)

 -------- Spline fit --------------------------- 

 x spline: x = [0.10558878272489601, 0.1305310750202113...
 y spline: y = [0.046372277538780926, 0.05201906296544364...

 ----------------------------------------------- 

```

Use `plot(fit.x,fit.y)` to plot the spline.












