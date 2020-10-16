#
# Perform multiple trials to find gloal minimum
#
function find_best_fit(f, ∇f, X, Y, np, options, lower, upper)
  local best_fit
  local fit
  func = (p) -> sq_residue(p, f, X, Y)
  grad! = (p,g) -> sq_residue_grad!(p,f,∇f,X,Y,g)
  best = +Inf
  nbest = 0
  ntrial = 1
  p0 = Vector{Float64}(undef,np)
  while nbest <= options.nbest && ntrial <= options.maxtrials
    ntrial += 1
    initP!(p0,options,lower,upper)
    if options.debug == false
      fit = SPGBox.spgbox!(p0, func, grad!, l=lower, u=upper) 
    else
      fit = SPGBox.spgbox!(p0, func, grad!, l=lower, u=upper, iprint=2) 
    end
    if fit.ierr == 0 
      if abs(fit.f - best) < options.besttol
        nbest = nbest + 1
        if fit.f < best
          best = fit.f
          best_fit = fit
        end
      elseif fit.f < best
        nbest = 1
        best = fit.f
        best_fit = fit
      end
    end
  end
  if nbest == 0
    error(" Could not obtain any successful fit, probably the data is not well posed.\n
            Further information can be obtained by adding Option(debug=true) to the input.")
  end
  return best_fit
end
