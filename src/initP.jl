#
# Random initial set of parameters
#
function initP!(p0::Vector{T}, options::Options, lower, upper) where {T}
    for i in eachindex(p0)
        pmin = max(T(options.p0_range[1]), lower[i])
        pmax = min(T(options.p0_range[2]), upper[i])
        p0[i] = pmin + rand() * (pmax - pmin)
    end
    return p0
end
