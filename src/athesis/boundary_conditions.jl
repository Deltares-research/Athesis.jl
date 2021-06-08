# Treatment of boundary conditions
# Using abstract loops on either CPU or GPU
using KernelAbstractions

function setBoundaries!(grid, state, bc)
    # For now we simply set 2 Dirichlet boundary conditions on the west and east boundaries
    # and we set 4 Neumann boundary conditions on the south, north, bottom and top boundaries
    # All boundaries are pressure boundaries for now

    deps = Event(device(grid.arch))

    # west and east
    event1 = launch!(
        grid, :yz, 
        westEastBoundary!, state.h, state.hⁿ⁺¹, grid.nx,
        bc.bcPressure[1], bc.bcPressure[2],
        dependencies=deps
        )

    # north and south
    event2 = launch!(
        grid, :xz, 
        northSouthBoundary!, state.h, state.hⁿ⁺¹, grid.ny,
        dependencies=deps
        )

    # bottom and top
    event3 = launch!(
        grid, :xy, 
        bottomTopBoundary!, state.h, state.hⁿ⁺¹, grid.nz,
        dependencies=deps
        )

    wait(device(grid.arch), event3)

    # Also set the four domain corners for plotting purposes
    # Simply copy one of the adjacent columns
    # These should by now be correctly set
    # hⁿ⁺¹[0,0,:] .= h[1,0,:]
    # hⁿ⁺¹[nx+1,0,:] .= h[nx,0,:]
    # hⁿ⁺¹[0,ny+1,:] .= h[0,ny,:]
    # hⁿ⁺¹[nx+1,ny+1,:] .= h[nx,ny+1,:]

end

@kernel function westEastBoundary!(h, hⁿ⁺¹, nx, bcWest, bcEast)
    j, k = @index(Global, NTuple)

    @inbounds begin
        h[0,j,k] = bcWest
        hⁿ⁺¹[0,j,k] = bcWest
        h[nx + 1,j,k] = bcEast
        hⁿ⁺¹[nx + 1,j,k] = bcEast
    end
end

@kernel function northSouthBoundary!(h, hⁿ⁺¹, ny)
    i, k = @index(Global, NTuple)

    @inbounds begin
        h[i,0,k] = h[i,1,k]
        hⁿ⁺¹[i,0,k] = h[i,1,k]
        h[i,ny + 1,k] = h[i,ny,k]
        hⁿ⁺¹[i,ny + 1,k] = h[i,ny,k]
    end
end

@kernel function bottomTopBoundary!(h, hⁿ⁺¹, nz)
    i, j = @index(Global, NTuple)

    @inbounds begin
    h[i,j,0] = h[i,j,1]
        hⁿ⁺¹[i,j,0] = h[i,j,1]
        h[i,j,nz + 1] = h[i,j,nz]
        hⁿ⁺¹[i,j,nz + 1] = h[i,j,nz]
    end
end

