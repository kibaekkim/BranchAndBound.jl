using Test

@testset "Unit test" begin
    include("unittest.jl")
end

@testset "milp example" begin
    include("../examples/milp.jl")
    @test tree.best_incumbent == -3.0
end

@testset "minlp example" begin
    include("../examples/minlp.jl")
    @test isapprox(tree.best_incumbent, -1)
end