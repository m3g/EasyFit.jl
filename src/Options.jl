@with_kw struct lower
  A = nothing
  B = nothing
  C = nothing
  a = nothing
  b = nothing
  c = nothing
  d = nothing
end

@with_kw struct upper
  A = nothing
  B = nothing
  C = nothing
  a = nothing
  b = nothing
  c = nothing
  d = nothing
end

@with_kw struct Options
  nbest :: Int = 5
  besttol :: Float64 = 1e-4
  maxtrials :: Int = 100
  fine :: Int = 100
  p0_range :: Vector{Float64} = [-1.,1.]
  lower :: lower = lower()
  upper :: upper = upper()
  debug :: Bool = false
end

Options(l :: lower) = Options(lower=l)
Options(u :: upper) = Options(upper=u)
Options(l :: lower, u :: upper) = Options(lower=l,upper=u)
Options(u :: upper, l :: lower) = Options(lower=l,upper=u)

export Options, lower, upper
