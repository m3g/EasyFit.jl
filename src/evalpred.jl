#
# Computes the prediction of a fit given the domain X
#
function evalpred(f,X,fit)
  ypred = Vector{Float64}(undef,length(X))
  for i in 1:length(X)
    ypred[i] = f(X[i],fit.x)
  end
  return ypred
end
