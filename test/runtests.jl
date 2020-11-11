using Athesis

using Test
using TimerOutputs

function testBulgeCPU()

    to = TimerOutput()

    input = getDefaultInput()
    bulge = input.Lx * input.Ly * input.const_recharge /
    (8.0 * input.K0 * input.boundary_pressure[1] * input.Δx * input.Δy)

    useGPU = false
    simulation = initSimulation(input, useGPU, to)
    runSimulation(simulation, to)

    ix = Int(ceil(simulation.grid.nx/2))
    iy = Int(ceil(simulation.grid.ny/2))
    iz = simulation.grid.nz
    head = simulation.state.h[ix,iy,iz]

    diff = abs(head - (input.boundary_pressure[1] + bulge))

    return diff < 0.0025
end

@testset "Athesis.jl" begin
    @testset "Validation tests" begin
        @test testBulgeCPU()
    end
end
