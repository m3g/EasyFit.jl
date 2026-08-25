```@meta
CurrentModule = EasyFit
```

# Bounds and fixed parameters

Lower and upper bounds can be set on the parameters of most functions using the
`l=lower(...)` and `u=upper(...)` keyword arguments. For example:

```@example bounds
using EasyFit, Random
Random.seed!(1)

x = sort(rand(10))
y = sort(rand(10))

fitlinear(x, y, l=lower(a=5.0), u=upper(a=10.0))
```

```@example bounds
y2 = @. 0.3 * exp(5x) + 0.7 * exp(-3x)
fitexp(x, y2, n=2, l=lower(a=[0.0, 0.0]), u=upper(a=[1.0, 1.0]))
```

Bounds on the intercepts or limiting values are not supported directly, but it is
possible to fix them to a constant value instead. For example:

```@example bounds
fitlinear(x, y, b=5.0)
```

```@example bounds
fitexp(x, y2, n=2, c=0.0)
```

The normalized exponential-decay fit, `fitexpdecay` (see [Normalized exponential decay](@ref)),
follows the same convention for its constant term `c`; its weights and decay rates
instead follow built-in constraints (``\sum_i a_i = 1`` and ``b_i > 0``) rather than
user-set bounds.
