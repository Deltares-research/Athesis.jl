using Test

function testBulgeCPU()
    testBulge(false)
end

function testBulgeGPU()
    	testBulge(true)
end

function testBulge(useGPU)
    to = TimerOutput()

    myFloat = Float64
    input = getDefaultInput(myFloat, "../Athesis.toml")

    # calculate bulging level from input parameters
    bulge = input.Lx * input.Lx * input.constRecharge /
    (8.0 * input.K0 * input.boundaryPressure[1] * input.Δx * input.Δy)

    useGPU = false
    simulation = initSimulation(input, useGPU, myFloat, to)
    runSimulation!(simulation, to)

    ix = Int(ceil(simulation.grid.nx/2))
    iy = Int(ceil(simulation.grid.ny/2))
    iz = simulation.grid.nz
    head = simulation.state.hⁿ⁺¹[ix,iy,iz]

    diff = abs(head - (input.boundaryPressure[1] + bulge))

    return diff < 1.0e-4
end

function testBulgeGPUvsCPU()

    to = TimerOutput()
    myFloat = Float64
    input = getDefaultInput(myFloat, "../Athesis.toml")

    # CPU
    simulationCPU = initSimulation(input, false, myFloat, to)
    convergedCPU, nIterCPU = runSimulation!(simulationCPU, to)

    # GPU
    simulationGPU = initSimulation(input, true, myFloat, to)
    convergedGPU, nIterGPU = runSimulation!(simulationGPU, to)

    # equal convergence
    @test convergedCPU == convergedGPU
    @test nIterCPU == nIterGPU

    # same result
    absDiffHeads = abs.(simulationCPU.state.hⁿ⁺¹ - simulationGPU.state.hⁿ⁺¹)
    maxDiff, = findmax(absDiffHeads)

    @test maxDiff ≈ 0.0 atol=1.0e-09

    return true
end


@testset "Validation" begin
    @info "Running validation tests on CPU..."
    @test testBulgeCPU()

    if has_cuda()
        @info "Running validation tests on GPU..."
        @test testBulgeGPU()
        @test testBulgeGPUvsCPU()
    end
end
