#
# Perform multiple trials to find gloal minimum
#
function find_best_fit(model, X, Y, np, options, lower, upper)
    local best_fit
    local fit
    best = +Inf
    nbest = 0
    ntrial = 1
    p0 = Vector{eltype(X)}(undef, np)
    while nbest <= options.nbest && ntrial <= options.maxtrials
        ntrial += 1
        try
            initP!(p0, options, lower, upper)
            fit = curve_fit(model, X, Y, p0, lower=lower, upper=upper)
            sum_residues = sum(abs2, fit.resid)
            fit_improve = sum_residues - best
            if fit_improve < options.besttol
                if fit_improve < 0
                    nbest = 1
                    best = sum_residues
                    best_fit = deepcopy(fit)
                else
                    nbest = nbest + 1
                end
            end
        catch msg
            if options.debug
                error("ERROR: $msg \n $fit")
            end
        end
    end
    if nbest == 0
        error("""
        Could not obtain any successful fit, probably the data is not well posed.
        Further information can be obtained by adding `options=Options(debug=true)` as kwarg.
        """)
    end
    return best_fit
end
