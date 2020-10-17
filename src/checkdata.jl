#
# Checks if the data is ok
#

function check_size(X,c)
  if (length(size(X)) > 2) || (length(size(X)) == 2 && size(X,2) != 1)
    error(" Only 1D arrays are accepted, and got $c with dimensions = ",size(X))
  end
  return vec(copy(X))
end

function checkdata(X :: Vectors, Y :: Vectors)
  if length(X) != length(Y)
    error(" Input x and y vectors must have the same length.")
  end
  X = check_size(X,"X")
  Y = check_size(Y,"Y")
  return X, Y
end
