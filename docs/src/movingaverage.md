```@meta
CurrentModule = EasyFit
```

# Moving averages

Use the `movavg` (or `movingaverage`) function. It computes the moving average of
`x[i]` in the range `i ± (n-1)/2` (if `n` is even, `n ← n + 1`):

```@example movingaverage
using EasyFit, Plots, Random
Random.seed!(1)

x = cumsum(randn(200))

ma = movavg(x, 21)
```

```@example movingaverage
plot(x, label="data", framestyle=:box, alpha=0.5)
plot!(ma.x, label="moving average", linewidth=2)
```
