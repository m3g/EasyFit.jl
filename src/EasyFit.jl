module EasyFit
  
  using Statistics
  using LsqFit
  using Interpolations
  using Parameters

  Numbers = Union{Int32,Int64,Float32,Float64}
  Vectors = Union{AbstractArray{Int32},
                  AbstractArray{Int64},
                  AbstractArray{Float32},
                  AbstractArray{Float64}}

  include("./setbounds.jl")
  include("./Options.jl")
  include("./checkdata.jl")
  include("./initP.jl")
  include("./finexy.jl")
  include("./pearson.jl")
  include("./find_best_fit.jl")

  include("./linear.jl")
  include("./quadratic.jl")
  include("./cubic.jl")
  include("./exponential.jl")
  include("./spline.jl")
  include("./movingaverage.jl")

end 
