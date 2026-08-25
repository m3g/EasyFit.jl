[![Tests](https://img.shields.io/badge/build-passing-green)](https://github.com/m3g/EasyFit.jl/actions)
[![codecov](https://codecov.io/gh/m3g/EasyFit.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/m3g/EasyFit.jl)
[![Aqua QA](https://JuliaTesting.github.io/Aqua.jl/dev/assets/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://m3g.github.io/EasyFit.jl/stable)
[![Docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://m3g.github.io/EasyFit.jl/dev)

# EasyFit

Easy interface for obtaining fits of 2D data.

The purpose of this package is to provide a very simple interface to obtain
some of the most common fits of 2D data: linear, quadratic, cubic, n-th degree
polynomial, exponential, normalized exponential-decay, and spline fits, plus
moving averages and density (continuous histogram) estimation.

On the background this interface uses the [LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl)
and [Interpolations](https://github.com/JuliaMath/Interpolations.jl) packages, which are
already quite easy to use. Additionally, EasyFit contains a simple globalization
heuristic, such that good non-linear fits are obtained often.

Our aim is to provide a package for quick fits without having to think about the code.

## Installation

```julia-repl
julia> ] add EasyFit

julia> using EasyFit
```

## Quick example

```julia-repl
julia> x = sort(rand(10)); y = sort(rand(10));

julia> fit = fitlinear(x, y)
------------------- Linear Fit -------------

Equation: y = ax + b

With: a = 1.158930569179642 ± 0.6538824813074927
      b = -0.1251714526967127 ± 0.3588742142656203

Correlation coefficient, R² = 0.9696101474036224
Average square residue = 0.004113279571449428

Predicted Y: ypred = [0.1044876257374221, 0.2072397615587609, ...]
residues = [0.08428396483020295, 0.05555828441380739, ...]

--------------------------------------------

julia> using Plots

julia> scatter(x, y)   # the original data

julia> plot!(fit.x, fit.y)   # the fit
```

## Documentation

The full documentation, with runnable examples and plots for every fit
(linear, quadratic, cubic, n-th degree polynomial, exponential, normalized
exponential-decay, splines, moving averages, density, bounds/fixed parameters,
and options) is available at:

**https://m3g.github.io/EasyFit.jl/stable**
