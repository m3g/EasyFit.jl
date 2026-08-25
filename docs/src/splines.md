```@meta
CurrentModule = EasyFit
```

# Splines

The fitting of splines requires the [Interpolations](https://github.com/JuliaMath/Interpolations.jl)
package to be loaded (this explicit requirement was introduced in version 0.6 of
`EasyFit`, and depends on `julia >= 1.9`).

Use the `fitspline` function:

```@example splines
using EasyFit, Interpolations, Plots, Random
Random.seed!(1)

x = sort(rand(10))
y = sort(rand(10))

fit = fitspline(x, y)
```

```@example splines
scatter(x, y, label="data", framestyle=:box)
plot!(fit.x, fit.y, label="spline", linewidth=2)
```
