#
# Perform multiple trials to find gloal minimum
#
function find_best_fit(model, X, Y, np, options)
  local best_fit
  best = +Inf
  nbest = 0
  ntrial = 1
  while nbest <= options.nbest && ntrial <= options.maxtrials
    ntrial += 1
    try 
      p0 = initP(np,options)
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
    end
  end
  return best_fit
end
