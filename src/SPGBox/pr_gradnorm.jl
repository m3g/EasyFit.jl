#
# Computes the norm of the projected gradient
# One method is defined for each possible scenario 
# of bounds (all, only lower, only upper, none). This is
# needed such that the lower and upper bounds vectors can 
# be of type Nothing if those bounds are not defined
#
function pr_gradnorm(x,g,l,u)
  if l == nothing && u == nothing 
    return pr_gradnorm_no_bounds(g)
  elseif l == nothing 
    pr_gradnorm_upper(x,g,u)
  elseif u == nothing 
    return pr_gradnorm_lower(x,g,l)
  else
    return pr_gradnorm_both_bounds(x,g,l,u)
  end
end
# Both bounds
function pr_gradnorm_both_bounds(x,g,l,u)
  gnorm = 0.
  for i in 1:length(x)
    z = max(l[i], min(u[i], x[i]-g[i])) - x[i]
    gnorm = max(gnorm, abs(z))
  end
  return gnorm
end
# Only lower
function pr_gradnorm_lower(x,g,l)
  gnorm = 0.
  for i in 1:length(x)
    z = max(l[i], x[i]-g[i]) - x[i]
    gnorm = max(gnorm, abs(z))
  end
  return gnorm
end
# Only upper
function pr_gradnorm_upper(x,g,u)
  gnorm = 0.
  for i in 1:length(x)
    z = min(u[i], x[i]-g[i]) - x[i]
    gnorm = max(gnorm, abs(z))
  end
  return gnorm
end
# No bounds
function pr_gradnorm_no_bounds(g)
  gnorm = 0.
  for i in 1:length(g)
    gnorm = max(gnorm, abs(g[i]))
  end
  return gnorm
end
