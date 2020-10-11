#
# Perform multiple trials to find gloal minimum
#
function find_best_fit(model, X, Y, np, options, lower, upper)
  local best_fit
  local fit
  best = +Inf
  nbest = 0
  ntrial = 1
  p0 = Vector{Float64}(undef,np)
  while nbest <= options.nbest && ntrial <= options.maxtrials
    ntrial += 1
    try 
      initP!(p0,options,lower,upper)
      fit = curve_fit(model, X, Y, p0, lower=lower, upper=upper)
      sum_residues = sum(fit.resid.^2)
      if abs(sum_residues - best) < options.besttol
        nbest = nbest + 1
        if sum_residues < best
          best = sum_residues
          best_fit = fit
        end
      elseif sum_residues < best
        nbest = 1
        best = sum_residues
        best_fit = fit
      end
    catch msg
      if options.debug
        error("ERROR: $msg \n $fit")
      end
    end
  end
  return best_fit
end
