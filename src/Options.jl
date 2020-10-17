@with_kw struct Options
  nbest :: Int = 5
  besttol :: Float64 = 1e-4
  maxtrials :: Int = 100
  fine :: Int = 100
  p0_range :: Vector{Float64} = [ -1e10 , 1e10 ]
  debug :: Bool = false
  spg_nitmax :: Int64 = 1000
  spg_nfevalmax :: Int64 = 10000
  spg_eps :: Float64 = 1.e-5
  spg_iprint :: Int64 = 0
end
export Options
