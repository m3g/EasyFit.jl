module EasyFit
  
  using Statistics
  using LsqFit
  using Interpolations
  using Parameters

  include("./LowerUpper.jl")
  include("./VarType.jl")
  include("./setbounds.jl")
  include("./Options.jl")
  include("./checkdata.jl")
  include("./initP.jl")
  include("./finexy.jl")
  include("./pearson.jl")
  include("./find_best_fit.jl")

  include("./FitMethods.jl")
  include("./linear.jl")
  include("./quadratic.jl")
  include("./cubic.jl")
  include("./exponential.jl")
  include("./spline.jl")
  include("./movingaverage.jl")

end 
