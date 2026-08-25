```@meta
CurrentModule = EasyFit
```

# Cubic fit

Use the `fitcubic` function:

```@example cubic
using EasyFit, Plots, Random
Random.seed!(1)

x = sort(rand(10))
y = x .^ 3 .+ 0.1 * rand(10)

fit = fitcubic(x, y)
```

```@example cubic
scatter(x, y, label="data", framestyle=:box)
plot!(fit.x, fit.y, label="fit", linewidth=2)
```
