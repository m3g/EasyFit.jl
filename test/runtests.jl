using TestItemRunner
@run_package_tests

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(EasyFit)
end