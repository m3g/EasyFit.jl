@with_kw struct Options
  nbest :: Int = 3
  besttol :: Float64 = 1e-4
  maxtrial :: Int = 100
  p0_range :: Vector{Float64} = [-1.,1.]
  fine :: Int = 100
end
export Options
