using Athesis

using Test
using TimerOutputs

function testBulgeCPU()

    to = TimerOutput()

    input = getDefaultInput()
    bulge = input.Lx * input.Lx * input.constRecharge /
    (8.0 * input.K0 * input.boundaryPressure[1] * input.Δx * input.Δy)

    useGPU = false
    simulation = initSimulation(input, useGPU, to)
    runSimulation!(simulation, to)

    ix = Int(ceil(simulation.grid.nx/2))
    iy = Int(ceil(simulation.grid.ny/2))
    iz = simulation.grid.nz
    head = simulation.state.hⁿ⁺¹[ix,iy,iz]

    diff = abs(head - (input.boundaryPressure[1] + bulge))

    return diff < 1.0e-4
end

@testset "Athesis.jl" begin
    @testset "Validation tests" begin
        @test testBulgeCPU()
    end
end
