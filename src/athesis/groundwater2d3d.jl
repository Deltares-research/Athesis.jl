# 3D groundwater

using Plots
using CuArrays

include("initialize.jl")
include("equations.jl")
include("external_forcing.jl")
include("postprocessing.jl")
include("model.jl")


# This is the present data storage:
# grid       = (ndim, nx, ny, Δx, Δy, x, y)
# state      = (h, u, v, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹)
# model      = externals
# parameters = (K, Δt, tend)
# externals  = source
# source     = (i_src, j_src, n_src, externals) or (i_src, j_src, k_src, n_src, externals)


function groundwater2d3d()

    # Initialize the model
    println("Groundwater 2D / 3D:")
    grid, model, state, parameters, time_data = model_initialize()

    # Unpack time parameters required for the time loop
    Δt       = time_data.Δt
    tend     = time_data.tend
    maxsteps = time_data.maxsteps
    time     = time_data.time

    # Set the plot title and layout
    #model_type = model.model_type
    #plot_title, plot_layout = set_plot_properties(model_type)
    #plotvars = prepare_plots(model_type, state)

    # Now start the time loop
    println("Starting time loop ...")
    time = Δt

    for n = 1:maxsteps

        # Add the sources
        set_sources!(time, model.source)
        pressure_equation!(grid, model, state, parameters, time_data)
        darcy_equation!(grid, model, state, parameters, time_data)

        # Set the new time
        time += Δt
        #println("Time step: ", n, ".      Time: ", time, " s.")

        if mod(n,1000)==0

            #println("Time step: ", n, ".      Time: ", time, " s.")
            #println("Plotting ...")
            plot_model(grid, state)

        end

    end
end

################################
# Execute the main function

@time groundwater2d3d()
################################
