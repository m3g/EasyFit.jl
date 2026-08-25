```@meta
CurrentModule = EasyFit
```

# N-th degree polynomial fit

Use the `fitndgr` function, passing the desired polynomial degree as the third argument:

```@example ndegree
using EasyFit, Plots, Random
Random.seed!(1)

x = sort(rand(10))
y = @. 1 + 2x + 3x^2 + 4x^3 + 6x^4

fit = fitndgr(x, y, 4)
```

The fitted coefficients `p[1], p[2], ..., p[n+1]` are stored in `fit.lscoeff`, from the
independent term to the highest-degree coefficient:

```@example ndegree
fit.lscoeff
```

```@example ndegree
scatter(x, y, label="data", framestyle=:box)
plot!(fit.x, fit.y, label="fit", linewidth=2)
```
