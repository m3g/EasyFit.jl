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
function fitspline(args...; kargs...)
    if isnothing(Base.get_extension(EasyFit, :SplineFitExt))
        error("Load first the `Interpolations` package to use the `fitspline` function.")
    end
    # The Interpolations extension is loaded but no method of fitspline matches
    # this call: raise a regular MethodError instead of a misleading message.
    throw(MethodError(fitspline, args))
end
@testitem "fitspline error" begin
    # Note: within the full test suite, other testitems already `using Interpolations`,
    # so the extension is loaded process-wide by the time this runs (package extensions,
    # once triggered, stay active for the rest of the session). The "Load first..."
    # message therefore can't be exercised here; it is covered by inspection/manual
    # testing in a fresh session instead. What we *can* verify in-process is the actual
    # bug fix: once the extension is loaded, a call that matches no method must raise a
    # regular (informative) MethodError instead of the misleading "Load first..." message.
    @test !isnothing(Base.get_extension(EasyFit, :SplineFitExt))
    @test_throws MethodError fitspline(1)
    @test_throws MethodError fitspline(1; x = 1)
    @test_throws MethodError fitspline(x = 1)
end

# fitexpdecay is defined in ext/ExpDecayFitExt.jl
export fitexpdecay
function fitexpdecay(args...; kargs...)
    if isnothing(Base.get_extension(EasyFit, :ExpDecayFitExt))
        error("Load first the `JuMP` and `Ipopt` packages to use the `fitexpdecay` function.")
    end
    # The JuMP/Ipopt extension is loaded but no method of fitexpdecay matches
    # this call: raise a regular MethodError instead of a misleading message.
    throw(MethodError(fitexpdecay, args))
end
@testitem "fitexpdecay error" begin
    # See the note in the "fitspline error" testitem above: the "Load first..."
    # message can't be exercised in-process here since other testitems already
    # load JuMP/Ipopt. What we verify here is the fix itself: once the extension
    # is loaded, a non-matching call raises a regular MethodError.
    @test !isnothing(Base.get_extension(EasyFit, :ExpDecayFitExt))
    @test_throws MethodError fitexpdecay(1)
    @test_throws MethodError fitexpdecay(1; n = 1)
end

end
