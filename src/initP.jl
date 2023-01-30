#
# Random initial set of parameters
#
function initP!(p0, options::Options, lower, upper)
    for i in eachindex(p0)
        pmin = max(options.p0_range[1], lower[i])
        pmax = min(options.p0_range[2], upper[i])
        p0[i] = pmin + rand() * (pmax - pmin)
    end
    return p0
end
