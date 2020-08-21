#
# Just returns a fine mesh for the fit
#
function finexy(X,fine,model,fit)
  Xmin = minimum(X)
  Xmax = maximum(X)
  x = collect(Xmin:(Xmax-Xmin)/fine:Xmax)
  y = model(x,fit.param)
  ypred = model(X,fit.param)
  return x, y, ypred
end
