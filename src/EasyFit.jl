module EasyFit
  
  using Statistics
  using LsqFit
  using Interpolations
  using Parameters

  Numbers = Union{Int64,Float64}
  Vectors = Union{AbstractArray{Int64},AbstractArray{Float64}}

  include("./Options.jl")
  include("./initP.jl")
  include("./finexy.jl")
  include("./pearson.jl")
  include("./find_best_fit.jl")

  include("./linear.jl")
  include("./quadratic.jl")
  include("./cubic.jl")
  include("./exponential.jl")
  include("./spline.jl")

end 
