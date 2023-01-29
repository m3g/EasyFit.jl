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
include("./fitlinear.jl")
include("./fitquadratic.jl")
include("./fitcubic.jl")
include("./fitpolynomial.jl")
include("./fitexponential.jl")
include("./fitspline.jl")
include("./movingaverage.jl")
include("./fitdensity.jl")

const model_catalogue = Dict(
  :linear => fitlinear,
  :quadratic => fitquadratic,
  :cubic => fitcubic,
  :polynomial => fitpolynomial,
  :exponential => fitexponential,
  :spline => fitspline,
  :movingaverage => movingaverage,
  :density => fitdensity,
)

end 
