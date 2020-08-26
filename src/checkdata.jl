#
# Checks if the data is ok
#
function checkdata(X :: Vectors, Y :: Vectors)
  if length(X) != length(Y)
    error(" Input x and y vectors must have the same length.")
  end
  return nothing
end
