#
# Checks if the data is ok
#

function check_size(X,c)
  if (length(size(X)) > 2) || (length(size(X)) == 2 && size(X,2) != 1)
    error(" Only 1D arrays are accepted, and got $c with dimensions = ",size(X))
  end
  return vec(copy(X))
end

function checkdata(X :: AbstractArray{<:Real}, Y :: AbstractArray{<:Real}, options :: Options)
  if length(X) != length(Y)
    error(" Input x and y vectors must have the same length.")
  end
  X = check_size(X,"X")
  Y = check_size(Y,"Y")
  # Set some reasonable ranges for the initial guesses
  options.p0_range[1] = minimum(Y) - 10*minimum(Y)
  options.p0_range[2] = maximum(Y) + 10*maximum(Y)
  return X, Y
end
