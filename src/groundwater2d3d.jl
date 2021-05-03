# 3D groundwater
using Plots
using CUDA
using TimerOutputs

using Athesis

function groundwater3d(isBenchmark = false, useGPU = false, myFloat = Float64)

    to = TimerOutput()

    @timeit to "total run time" begin

        # Initialize the model
        println("Running 3D groundwater model:")

        defaultInput = getDefaultInput(myFloat)
        simulation = initSimulation(defaultInput, useGPU, myFloat, to)

        runSimulation!(simulation, to)

        plotSimulation(simulation, to)

    end

    print_timer(to)

end
