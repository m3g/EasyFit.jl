module EasyFit

using TestItems
using Statistics
using LsqFit
using Parameters
using Unitful: ustrip, oneunit, unit

# supertype for all fits, to help on dispatch of common methods
abstract type Fit{T<:AbstractFloat} end

include("./LowerUpper.jl")
include("./VarType.jl")
include("./setbounds.jl")
include("./Options.jl")
include("./checkdata.jl")
include("./initP.jl")
include("./finexy.jl")
include("./R2.jl")
include("./find_best_fit.jl")

include("./FitMethods.jl")
include("./fitlinear.jl")
include("./fitquadratic.jl")
include("./fitcubic.jl")
include("./fitndgr.jl")
include("./fitexponential.jl")
include("./movingaverage.jl")
include("./fitdensity.jl")

# fitspline is defined in ext/SplineFitExt.jl
export fitspline
fitspline(args...; kargs...) = error("Load first the `Interpolations` package to use the `fitspline` function.")
@testitem "fitspline error" begin
    @test_throws "Load first the `Interpolations` package to use the `fitspline` function." fitspline(1)
    @test_throws "Load first the `Interpolations` package to use the `fitspline` function." fitspline(1; x = 1)
    @test_throws "Load first the `Interpolations` package to use the `fitspline` function." fitspline(x = 1)
end

# fitexpdecay is defined in ext/ExpDecayFitExt.jl
export fitexpdecay
fitexpdecay(args...; kargs...) = error("Load first the `JuMP` and `Ipopt` packages to use the `fitexpdecay` function.")
@testitem "fitexpdecay error" begin
    @test_throws "Load first the `JuMP` and `Ipopt` packages to use the `fitexpdecay` function." fitexpdecay(1)
    @test_throws "Load first the `JuMP` and `Ipopt` packages to use the `fitexpdecay` function." fitexpdecay(1; n = 1)
end

end
