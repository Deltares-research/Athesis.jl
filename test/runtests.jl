using Athesis

using Test

include("../src/groundwater2d3d.jl")

@testset "Athesis.jl" begin
    @testset "Unit tests" begin
        @test 2 + 2 == 4
    end

    @testset "Performance tests" begin
        groundwater3d(true,false)
        @test true
    end
end
