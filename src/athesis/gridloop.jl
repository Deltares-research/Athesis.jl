using CUDA

# This is the present data storage:
# grid       = (nx, ny, nz, Δx, Δy, Δz, x, y, z)
# state      = (h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
# model      = externals
# parameters = (K, Δt, tend)
# externals  = source
# source     = (i_src, j_src, k_src, n_src, externals)


function launch!(grid, dims, kernel!, args...; dependencies=nothing)

    if dims == :xyz
        workgroup = (8,8,4)
        worksize = (grid.nx, grid.ny, grid.nz)
    elseif dims == :xy
        workgroup = (16,16)
        worksize = (grid.nx, grid.ny)
    elseif dims == :yz
        workgroup = (16,16)
        worksize = (grid.ny, grid.nz)
    elseif dims == :xz
        workgroup = (16,16)
        worksize = (grid.nx, grid.nz)
    else
        throw(ArgumentError("Error launching kernel, dims: $dims"))
    end



    @debug "Launching kernel $kernel! with worksize $worksize"

    _kernel! = kernel!(device(grid.arch), workgroup, worksize)
    event = _kernel!(args..., dependencies=dependencies)

    return event
end
