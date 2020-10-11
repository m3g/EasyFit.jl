#
# Function to set bounds with optional user input 
#

function lower_bound(value,varname,vartype) 
  value == nothing && return -Inf
  value isa vartype || error("Lower bound of $varname must be of type $vartype")
  return value
end

function upper_bound(value,varname,vartype) 
  value == nothing && return +Inf
  value isa vartype || error("Upper bound of $varname must be of type $vartype")
  return value
end
