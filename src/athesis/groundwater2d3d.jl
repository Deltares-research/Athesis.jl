# 3D groundwater

using Plots
using CuArrays

include("initialize.jl")
include("equations.jl")
include("external_forcing.jl")
include("postprocessing.jl")
include("model.jl")


# This is the present data storage:
# grid       = (nx, ny, nz, Δx, Δy, Δz, x, y, z)
# state      = (h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
# model      = externals, recharge
# parameters = (K, Δt, tend)
# externals  = source
# source     = (i_src, j_src, k_src, n_src, externals)


function groundwater3d()

    # Initialize the model
    println("Running 3D groundwater model:")
    grid, model, state, parameters, time_data = model_initialize()

    # Unpack time parameters required for the time loop
    Δt       = time_data.Δt
    tend     = time_data.tend
    maxsteps = time_data.maxsteps
    time     = time_data.time

    # Now start the time loop
    println("Starting time loop ...")

    # Initialize the present time
    time = 0.0

    for n = 1:maxsteps

        # Add the sources
        set_sources!(time, model.source)
        set_recharge!(time,model.recharge)
        pressure_equation!(grid, model, state, parameters, time_data)
        darcy_equation!(grid, model, state, parameters, time_data)

        # Set the new time
        time += Δt
        #println("Time step: ", n, ".      Time: ", time, " s.")

        if mod(n,50)==0

            #println("Time step: ", n, ".      Time: ", time, " s.")
            #println("Plotting ...")
            plot_model(grid, state)

        end

    end
end

################################
# Execute the main function

@time groundwater3d()
################################
