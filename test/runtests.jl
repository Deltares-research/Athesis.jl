using Athesis

using Test

include("../src/groundwater2d3d.jl")

@testset "Athesis.jl" begin
    @testset "Unit tests" begin
        @test 2 + 2 == 4
    end

    @testset "Performance tests" begin
        #groundwater3d(true,false)

        opbolling = 80.0*80.0*5.0e-4/(8.0*10.0)

        #println()
        #println(opbolling)
        #println()
        @test true
    end
end
