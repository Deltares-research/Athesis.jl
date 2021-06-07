# equations.jl

# This file contains all implemented equations solved in the model

# Currently implemented:
# - Pressure equation
# - Darcy equation (velocity)

# This is the present data storage:
# grid       = (nx, ny, nz, Δx, Δy, Δz, x, y, z)
# state      = (h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
# model      = externals, recharge
# parameters = (K, Δt, tend)
# externals  = source
# source     = (i_src, j_src, k_src, n_src, externals)


function pressureEquation!(grid, model, state, parameters, timeData)

    source = model.source.externalSource
    rchFlux = model.recharge.rechargeFlux

    # Solve for the pressure/head
    deps = Event(device(grid.arch))

    toplayerNr = grid.nz
    event1 = launch!(
                    grid,               # grid
                    :xy,                # coverage
                    addToLayer!,        # kernel                    
                    source,             # kernel args
                    rchFlux,
                    toplayerNr,
                    dependencies=deps   # handle for sync.
                    )
                    
    event2 = launch!(
                    grid,                # grid
                    :xyz,                # coverage
                    pressureKernel!,  # kernel
                    source,              # kernel args...
                    state.h,
                    state.hⁿ⁺¹,
                    grid.Δx, 
                    grid.Δy, 
                    grid.Δz, 
                    timeData.Δt, 
                    parameters.K,
                    parameters.specificStorage, 
                    dependencies=event1  # handle for sync.
                    )

    wait(device(grid.arch), event2)

end

function darcyEquation!(grid, model, state, parameters, timeData)
    # Compute the velocities from the darcy equation using the pressure
    deps = Event(device(grid.arch))
    event = launch!(
                    grid,                # grid
                    :xyz,                # coverage
                    darcyKernel!,     # kernel
                    state.uⁿ⁺¹,          # kernel args...
                    state.vⁿ⁺¹,
                    state.wⁿ⁺¹,
                    state.h,
                    grid.Δx, 
                    grid.Δy, 
                    grid.Δz, 
                    parameters.K,
                    dependencies=deps # event handle for sync.
                    )
    wait(device(grid.arch), event)
end
