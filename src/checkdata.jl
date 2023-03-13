#
# Checks if the data is ok
#
function check_size_and_type(X, c, data_type)
    if (length(size(X)) > 2) || (length(size(X)) == 2 && size(X, 2) != 1)
        throw(ArgumentError("Only 1D arrays are accepted, and got $c with dimensions = $size(X)"))
    end
    return vec(data_type.(X))
end

function checkdata(X::AbstractArray{T1}, Y::AbstractArray{T2}, options::Options) where {T1<:Real, T2<:Real}
    if length(X) != length(Y)
        throw(ArgumentError("Input x and y vectors must have the same length."))
    end
    data_type = promote_type(Float32,T1,T2)
    X = check_size_and_type(X, "X", data_type)
    Y = check_size_and_type(Y, "Y", data_type)
    # Set some reasonable ranges for the initial guesses
    if options.p0_range[1] ≈ 0.0 && options.p0_range[2] ≈ 0.0
        range = maximum(Y) - minimum(Y)
        options.p0_range[1] = -100 * range
        options.p0_range[2] = +100 * range
    end
    return X, Y, data_type
end
