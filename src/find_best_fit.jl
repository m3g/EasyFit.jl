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
  auxvecs = SPGBox.VAux(np)
  while nbest <= options.nbest && ntrial <= options.maxtrials
    ntrial += 1
    initP!(p0,options,lower,upper)
    if options.debug
      println(" Initial parameters: ", p0 )
    end
    fit = SPGBox.spgbox!(p0, func, grad!, l=lower, u=upper,
                         vaux = auxvecs,
                         iprint = options.spg_iprint,
                         eps = options.spg_eps,
                         nitmax = options.spg_nitmax, 
                         nfevalmax = options.spg_nfevalmax) 
    if options.debug
      println(" Trial number = ", ntrial)
      println(fit)
    end
    if fit.ierr == 0 
      if abs(fit.f - best) < options.besttol
        nbest = nbest + 1
        if fit.f < best
          best = fit.f
          best_fit = deepcopy(fit)
        end
      elseif fit.f < best
        nbest = 1
        best = fit.f
        best_fit = deepcopy(fit)
      end
    end
  end
  if nbest == 0
error(" Could not obtain any successful fit, probably the data is not well posed.
        Further information can be obtained by adding options=Option(debug=true) 
        to the input.")
  end
  return best_fit
end
