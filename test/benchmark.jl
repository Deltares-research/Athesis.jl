using Test

@testset "Benchmark" begin
    @info "Running benchmark tests..."

    include("../benchmark/benchmark_hooghoudt.jl")
end
