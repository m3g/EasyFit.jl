#
# Returns the R^2 correlation coefficient
#
function R2(X, Y, model, fit)
    ypred = similar(Y)
    ypred = model(X, fit.param)
    R = Statistics.cor(Y, ypred)
    return R^2
end
