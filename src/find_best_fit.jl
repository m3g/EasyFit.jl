#
# Perform multiple trials to find gloal minimum
#
function find_best_fit(model, X, Y, np, options)
  local best_fit
  local fit
  best = +Inf
  nbest = 0
  ntrial = 1
  p0 = initP(np,options)
  while nbest <= options.nbest && ntrial <= options.maxtrials
    ntrial += 1
    try 
      nextP!(p0,options)
      fit = curve_fit(model, X, Y, p0)
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
