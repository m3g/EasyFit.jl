```@meta
CurrentModule = EasyFit
```

# Normalized exponential decay

`fitexpdecay` fits a **normalized** multiple-exponential decay:

```math
y(t) = \sum_{i=1}^n a_i \, e^{-t/b_i} + c
```

subject to the constraints ``\sum_i a_i + c = 1`` (so that ``y(0) = 1`` always, by
construction), ``b_i > 0`` for every decay rate, and (when `c` is fitted freely)
``c \geq 0``. This is the common form used, for example, to describe fluorescence or
other time-resolved decays where the amplitudes are fractional populations, and the
constant ``c`` is the long-time baseline fraction — together they make up the whole
population, which is normalized to one at ``t=0``.

Because the amplitudes and time constants are optimized subject to nonlinear
constraints, this fit is implemented as a
[package extension](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages)
that requires [JuMP](https://github.com/jump-dev/JuMP.jl) and
[Ipopt](https://github.com/jump-dev/Ipopt.jl) to be loaded:

```@example expdecay
using EasyFit, JuMP, Ipopt, Plots, Random
Random.seed!(1)
```

## Basic usage

`fitexpdecay` accepts a single vector of values. By default, the time associated to
each data point is simply its index step, i.e. the first data point is assumed to
correspond to `t = 0`:

```@example expdecay
t = 0:0.1:8
y = @. 0.7 * exp(-t / 0.6) + 0.3 * exp(-t / 4) + 0.01 * randn()

fit = fitexpdecay(y; n=2, t=collect(t))
```

The fitted weights and constant always add up to one (so `fit(0) == 1`), and the
decay rates are always positive:

```@example expdecay
sum(fit.a) + fit.c, all(fit.b .> 0)
```

```@example expdecay
scatter(t, y, label="data", framestyle=:box, markersize=3, markerstrokewidth=0)
plot!(fit.x, fit.y, label="fit", linewidth=2)
```

## Index-based time and `OffsetArray`s

If no time vector is given, `fitexpdecay` uses the position of each data point
(starting at zero) as its time, i.e. consecutive samples are one time unit apart.
This is convenient for data naturally indexed by an `OffsetArray` starting at `0`:

```@example expdecay
using OffsetArrays

yindex = OffsetArray([0.7 * exp(-i / 6) + 0.3 * exp(-i / 40) for i in 0:79], 0:79)
fit_default_times = fitexpdecay(yindex; n=2)

fit_default_times.a, fit_default_times.b
```

which recovers the same decay rates (`6` and `40`) that were used to generate the
data on its index scale, without having to build a separate time vector.

## Fixing the constant term

Because ``\sum_i a_i + c = 1`` always holds, `y(0) = 1` regardless of how `c` is
chosen — fixing `c` only fixes the split between the long-time baseline fraction
(`c`) and the fraction explained by decay (`sum(a) = 1 - c`), it does not change the
overall normalization. The independent constant `c` is fitted freely by default,
subject to ``c \geq 0``. It can optionally be fixed to a user-provided value with the
`c` keyword — a fixed value is *not* required to be non-negative:

```@example expdecay
y_baseline = @. 0.3 + 0.7 * exp(-t / 2) + 0.01 * randn()
fit_fixed_c = fitexpdecay(y_baseline; n=1, t=collect(t), c=0.3)

fit_fixed_c.c, sum(fit_fixed_c.a)
```
