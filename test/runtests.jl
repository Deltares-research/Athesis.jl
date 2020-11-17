using Test
using TimerOutputs
using CUDA

using Athesis

@testset "Athesis tests" begin
    include("unit.jl")
    include("validation.jl")
    include("benchmark.jl")
end
