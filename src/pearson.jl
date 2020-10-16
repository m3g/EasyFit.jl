#
# Returns the pearson coefficient of the fit
#
function pearson(X,Y,f,fit)
  ypred = evalpred(X,fit)
  R = Statistics.cor(Y,ypred)
  return R
end
