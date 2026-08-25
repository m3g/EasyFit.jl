```@meta
CurrentModule = EasyFit
```

# Quadratic fit

Use the `fitquad` (or `fitquadratic`) function:

```@example quadratic
using EasyFit, Plots, Random
Random.seed!(1)

x = sort(rand(10))
y = x .^ 2 .+ 0.1 * rand(10)

fit = fitquad(x, y)
```

```@example quadratic
scatter(x, y, label="data", framestyle=:box)
plot!(fit.x, fit.y, label="fit", linewidth=2)
```
