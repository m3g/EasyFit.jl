#
# Checks if the data is ok
#
function check_size(X, c)
    if (length(size(X)) > 2) || (length(size(X)) == 2 && size(X, 2) != 1)
        throw(ArgumentError("Only 1D arrays are accepted, and got $c with dimensions = $size(X)"))
    end
    return vec(copy(X))
end

function checkdata(X::AbstractArray{<:Real}, Y::AbstractArray{<:Real}, options::Options)
    if length(X) != length(Y)
        throw(ArgumentError("Input x and y vectors must have the same length."))
    end
    X = check_size(X, "X")
    Y = check_size(Y, "Y")
    # Set some reasonable ranges for the initial guesses
    if options.p0_range[1] ≈ 0.0 && options.p0_range[2] ≈ 0.0
        range = maximum(Y) - minimum(Y)
        options.p0_range[1] = -100 * range
        options.p0_range[2] = +100 * range
    end
    return X, Y
end
