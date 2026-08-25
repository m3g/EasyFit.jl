```@meta
CurrentModule = EasyFit
```

# Linear fit

To perform a linear fitting, use the `fitlinear` function:

```@example linear
using EasyFit, Plots, Random
Random.seed!(1)

x = sort(rand(10))
y = sort(rand(10)) # some data

fit = fitlinear(x, y)
```

The `fit` data structure which comes out of `fitlinear` contains the output data with
the same names as shown in the above output:

```@example linear
fit.a, fit.sd_a
```

```@example linear
fit.b, fit.sd_b
```

```@example linear
fit.R2
```

The `fit.x` and `fit.y` vectors can be used for plotting the results, together with the
original data:

```@example linear
scatter(x, y, label="data", framestyle=:box)
plot!(fit.x, fit.y, label="fit", linewidth=2)
```
