#
# Just returns a fine mesh for the fit
#
function finexy(f,X,fine,fit)
  Xmin = minimum(X)
  Xmax = maximum(X)
  x = collect(Xmin:(Xmax-Xmin)/fine:Xmax)
  y = evalpred(f,x,fit)
  return x, y
end
