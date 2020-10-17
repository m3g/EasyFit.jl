#
# Returns the pearson coefficient of the fit
#
function pearson(X,Y,model,fit)
  ypred = similar(Y) 
  ypred = model(X,fit.param)
  R = Statistics.cor(Y,ypred)
  return R
end
