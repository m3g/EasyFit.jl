module EasyFit
  
  using Statistics
  using Interpolations
  using Parameters

  include("./SPGBox/SPGBox.jl")
  include("./sq_residue.jl")

  include("./evalpred.jl")
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
  include("./fitlinear.jl")
  include("./fitquadratic.jl")
  include("./fitcubic.jl")
  include("./fitexponential.jl")
  include("./fitspline.jl")
  include("./movingaverage.jl")

end 
