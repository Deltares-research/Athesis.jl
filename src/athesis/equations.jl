# equations.jl

# This file contains all implemented equations solved in the model

# Currently implemented:
# - Pressure equation
# - Darcy equation (velocity)

# This is the present data storage:
# grid       = (ndim, nx, ny, Δx, Δy, x, y) or (ndim, nx, ny, nz, Δx, Δy, Δz, x, y, z)
# state      = (h, u, v, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹) or (h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
# model      = externals
# parameters = (K, Δt, tend)
# externals  = source
# source     = (i_src, j_src, n_src, externals) or (i_src, j_src, k_src, n_src, externals)

include("gridloop.jl")
include("kernels.jl")


function pressure_equation!(grid, model, state, parameters, time_data)

    # 3D implementation
    #println("   Solving the pressure equation ...")

    # Unpack only source to be able to dispatch on type (Array or CuArray)
    source = model.source.external_source

    # Solve for the pressure/head
    gridloop(pressure_kernel, source, state, grid, parameters, time_data)

    # Update the old to the new solution
    state.h = copy(state.hⁿ⁺¹)
end



function darcy_equation!(grid, model, state, parameters, time_data)

    # Compute the velocities from the darcy equation using the pressure
    # For now cell-centered and collocated
    #println("   Solving the Darcy equation ...")

    # Unpack only source to be able to dispatch on type (Array or CuArray)
    source = model.source.external_source

    gridloop(darcy_kernel, source, state, grid, parameters, time_data)

    # Update the old to the new solution
    state.u = copy(state.uⁿ⁺¹)
    state.v = copy(state.vⁿ⁺¹)
end
