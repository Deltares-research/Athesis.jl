# This script contains the plotting of timer output
using Plots

function plotTimerOutput(archs, dims, to::TimerOutput, modelName)

    x = fill(0, length(dims))
    y = fill(0.0, length(dims))

    # Build a list of architectures in a single string, to include in the title
    archList = ""
    for arch in archs
        if arch != archs[end]
            archList *= arch*" vs. "
        else
            archList *= arch*" "
        end
    end

    # Make the plot title
    plotTitle = "Performance comparison for " * archList * "for "* modelName * " model"

    # Start the empty plot with title
    p = plot(title=plotTitle)

    # Create placeholder for the reference time (used for computing the speed-up)
    refTime = 0.0

    # For now allow up to 5 different architectures.
    # When more are needed, expand this list with more markers.
    mark = [(:circle), (:square), (:cross), (:pentagon), (:hexagon)]

    # Now loop over all architecture and grid resolutions to build up the plot
    archIdx = 0
    for arch in archs
        archIdx += 1
        dimIdx = 0
        for dim in dims
            dimIdx += 1
            dimstr = string(dim)
            time = TimerOutputs.time(to[arch][dimstr])

            # Determine total number of cells for this resolution
            numcells = 1
            for idim in dim
                numcells *=idim
            end
            x[dimIdx] = numcells
            # Time was given in nano-seconds: convert to seconds
            y[dimIdx] = Float64(time)/1.0e9
        end
        if arch == archs[1]
            refTime = copy(y)
        end
        p = plot!(x,y, marker=mark[archIdx], size=(800,800), xlabel="# of cells", ylabel="Time (s)", label=arch, legend=(0.2, 0.9), grid=true)
    end

    @show SpeedUp = refTime./y
    p = plot!(x,SpeedUp, marker=:diamond, size=(800,800), xlabel="# of cells", ylabel="Time (s)", label="Architectural speed up", legend=(0.2, 0.9), grid=true)

    #display(p)
    figureTitle = "performance_plot_" * modelName * ".png"
    #savefig(p, figureTitle)


end
