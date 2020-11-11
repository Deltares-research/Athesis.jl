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

    #println("   Solving the pressure equation ...")

    # Unpack only source to be able to dispatch on type (Array or CuArray)
    source   = model.source.externalSource

    # Add the model recharge to top layer
    source[:,:,grid.nz] .+= model.recharge.rechargeFlux

    # Solve for the pressure/head
    gridloop!(pressureKernel!, source, state, grid, parameters, timeData)
end



function darcyEquation!(grid, model, state, parameters, timeData)

    # Compute the velocities from the darcy equation using the pressure
    # For now cell-centered and collocated
    #println("   Solving the Darcy equation ...")

    # Unpack only source to be able to dispatch on type (Array or CuArray)
    source = model.source.externalSource

    gridloop!(darcyKernel!, source, state, grid, parameters, timeData)
end
