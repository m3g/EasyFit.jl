import CompatHelperLocal as CHL
CHL.@check()
using EasyFit
using Random
using Statistics
using Test

Random.seed!(0)

using Logging
is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
if is_logging(stderr)
    global_logger(NullLogger())
end

@testset "EasyFit.jl" begin

    x = sort(rand(10))
    y = rand(10)

    for (model, fun) in model_catalogue
        f = fun(x, y)
        @testset "Estimates" begin
            @test all(y - f.ypred .== f.residues)
            ss_res = sum(f.residues.^2)
            ss_tot = sum((y .- mean(y)).^2)
            @test isapprox(f.R^2, 1 - (ss_res/ss_tot), atol=1e-5)
            @test all(f.ypred == f.(x))
        end
    end

end
