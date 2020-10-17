@with_kw struct Options
  nbest :: Int = 5
  besttol :: Float64 = 1e-4
  maxtrials :: Int = 100
  fine :: Int = 100
  p0_range :: Vector{Float64} = [-1.,1.]
  rand_std :: Float64 = 0.5
  debug :: Bool = false
end
export Options
