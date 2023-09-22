module EasyFit

using TestItems
using Statistics
using LsqFit
using Parameters

# supertype for all fits, to help on dispatch of common methods
abstract type Fit{T<:AbstractFloat} end

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
include("./movingaverage.jl")
include("./fitdensity.jl")

# fitspline is defined in ext/SplineFitExt.jl
export fitspline
fitspline(args...) = error("Load first the `Interpolations` package to use the `fitspline` function.")
@testitem "fitspline error" begin
    @test_throws "Load first the `Interpolations` package to use the `fitspline` function." fitspline()
end

end
