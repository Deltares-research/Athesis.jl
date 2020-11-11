# 3D groundwater
using Plots
using CUDA
using TimerOutputs

using Athesis

function groundwater3d(isBenchmark = false, useGPU = false)

    to = TimerOutput()

    @timeit to "total run time" begin

        # Initialize the model
        println("Running 3D groundwater model:")

        @timeit to "initialization" CUDA.@sync begin
            defaultInput = getDefaultInput()
            simulation = initSimulation(defaultInput, useGPU, to)
        end

        @timeit to "run"  CUDA.@sync begin
            runSimulation!(simulation, to)
        end

        @timeit to "plot result" CUDA.@sync begin
            plotSimulation(simulation)
        end
    end

    print_timer(to)

end
