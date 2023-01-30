#
# Just returns a fine mesh for the fit
#
function finexy(X, fine, model, fit)
    Xmin, Xmax = extrema(X)
    x = collect(range(Xmin, Xmax, length=fine))
    y = model(x, fit.param)
    ypred = model(X, fit.param)
    return x, y, ypred
end
