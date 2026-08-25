```@meta
CurrentModule = EasyFit
```

# Exponential fits

Use the `fitexp` (or `fitexponential`) function for a single exponential:

```@example exponential
using EasyFit, Plots, Random
Random.seed!(1)

x = sort(rand(10))
y = @. 0.3 * exp(5x) + 0.1 * rand()

fit = fitexp(x, y)
```

```@example exponential
scatter(x, y, label="data", framestyle=:box)
plot!(fit.x, fit.y, label="mono-exponential fit", linewidth=2)
```

Add `n=N` for a sum of `N` exponentials:

```@example exponential
y2 = @. 0.3 * exp(5x) + 0.7 * exp(-3x) + 0.05 * rand()

fit2 = fitexp(x, y2, n=2)
```

```@example exponential
scatter(x, y2, label="data", framestyle=:box)
plot!(fit2.x, fit2.y, label="bi-exponential fit", linewidth=2)
```

The intercept `c` is fitted freely by default; it can be fixed to a constant value
with the `c` keyword, and lower/upper bounds can be set as described in
[Bounds and fixed parameters](@ref):

```@example exponential
fitexp(x, y2, n=2, c=0.0)
```
