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

include("gridloop.jl")
include("kernels.jl")


function pressure_equation!(grid, model, state, parameters, time_data)

    #println("   Solving the pressure equation ...")

    # Unpack only source to be able to dispatch on type (Array or CuArray)
    source   = model.source.external_source

    # Add the model recharge to top layer
    source[:,:,grid.nz] .+= model.recharge.recharge_flux

    # Solve for the pressure/head
    gridloop!(pressure_kernel!, source, state, grid, parameters, time_data)
end



function darcy_equation!(grid, model, state, parameters, time_data)

    # Compute the velocities from the darcy equation using the pressure
    # For now cell-centered and collocated
    #println("   Solving the Darcy equation ...")

    # Unpack only source to be able to dispatch on type (Array or CuArray)
    source = model.source.external_source

    gridloop!(darcy_kernel!, source, state, grid, parameters, time_data)
end
