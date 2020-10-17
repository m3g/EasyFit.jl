#
# Perturbed set of parameters for next trial
#
function nextP!(p,options)
  @. p = p + (options.p0_range[2]-options.p0_range[1])*(randn()*options.rand_std)
end

