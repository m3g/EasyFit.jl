#
# Random initial set of parameters
#
function initP(n :: Int64, options :: Options)
  p0 = Vector{Float64}(undef,n)
  @. p0 = options.p0_range[1] + (options.p0_range[2]-options.p0_range[1])*rand()
  return p0
end
