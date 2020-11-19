# this script runs the Hooghoudt model benchmark
using Athesis

using TimerOutputs

include("../benchmark/benchmark_utils.jl")

to = TimerOutput()

archs = ["CPU"]
@withCUDA archs = ["CPU", "GPU"]

dims = [(16, 16, 16), (32, 32, 16), (64, 64, 16), (128, 128, 16)] #, (256, 256, 16), (512, 512, 16)]

for arch in archs
    @synctimeit to arch begin
        for dim in dims
            input = getDefaultInput()
            input.nx = dim[1]
            input.ny = dim[2]
            input.nz = dim[3]

            # enforce precompilation
            @info "starting precompilation..."
            preSim = initSimulation(input, arch=="GPU")
            doTimestep!(1, preSim)

            # run timing
            label = string(dim)
            @info "run benchmark: " * label
            @synctimeit to label begin
                sim = initSimulation(input, arch == "GPU")

                nmax = 200
                for i in 1:nmax
                    hasCvg, = doTimestep!(i, sim)
                    if hasCvg
                        break
                    end
                end

                @test true
            end
        end
    end
end

print_timer(to)

modelName = "Hooghoudt"
plotTimerOutput(archs, dims, to, modelName)
