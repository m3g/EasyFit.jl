#
# Structure to store auxiliary arrays
#
struct Aux
  g :: Vector{Float64}
  xn :: Vector{Float64}
  gn :: Vector{Float64}
  fprev :: Vector{Float64}
end
Aux(n,m) = Aux( Vector{Float64}(undef,n), 
                Vector{Float64}(undef,n), 
                Vector{Float64}(undef,n), 
                Vector{Float64}(undef,m) )
# The default value for m
Aux(n) = Aux(n,10)
