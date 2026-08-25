```@meta
CurrentModule = EasyFit
```

# Options

It is possible to pass an optional set of parameters to the fitting functions with the
`options` keyword, using an `Options` struct. For example:

```@example options
using EasyFit, Random
Random.seed!(1)

x = sort(rand(10))
y = @. 0.3 * exp(5x)

fitexp(x, y, options=Options(maxtrials=1000))
```

Available options:

| Keyword | Type | Default value | Meaning |
|:-------:|:----:|:-------------:|:--------|
| `fine`  | `Int`| 100           | Number of points of fit to smooth plot. |
| `p0_range`  | `Vector{Float64}`  | `[-100*(maximum(Y)-minimum(Y)), 100*(maximum(Y)-minimum(Y))]`  | Range of generation of initial random parameters. |
| `nbest` | `Int`| 5  | Number of repetitions of best solution in global search. |
| `besttol` | `Float64`| 1e-4  | Similarity of the sum of residues of two solutions such that they are considered the same. |
| `maxtrials`  | `Int`| 100  | Maximum number of trials in global search. |
| `debug` | `Bool` | false | Prints errors of failed fits. |

The `nbest`, `besttol`, and `maxtrials` options are also used by the multistart search
performed by `fitexpdecay` (see [Normalized exponential decay](@ref)).
