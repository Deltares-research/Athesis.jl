# Treatment of boundary conditions
# Using abstract loops on either CPU or GPU

function set_boundaries!(grid, state, bc)
    # For now we simply set 2 Dirichlet boundary conditions on the west and east boundaries
    # and we set 4 Neumann boundary conditions on the south, north, bottom and top boundaries
    # All boundaries are pressure boundaries for now

    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    boundaryloop!(westBoundary!, state.h, state.hⁿ⁺¹, ny, nz, nx+1, bc.bc_pressure[1])

    boundaryloop!(eastBoundary!, state.h, state.hⁿ⁺¹, ny, nz, nx+1, bc.bc_pressure[2])

    boundaryloop!(southBoundary!, state.h, state.hⁿ⁺¹, nx, nz, ny+1, 0.0)

    boundaryloop!(northBoundary!, state.h, state.hⁿ⁺¹, nx, nz, ny+1, 0.0)

    boundaryloop!(bottomBoundary!, state.h, state.hⁿ⁺¹, nx, ny, nz+1, 0.0)

    boundaryloop!(topBoundary!, state.h, state.hⁿ⁺¹, nx, ny, nz+1, 0.0)

    # Also set the four domain corners for plotting purposes
    # Simply copy one of the adjacent columns
    # These should by now be correctly set
    # hⁿ⁺¹[0,0,:] .= h[1,0,:]
    # hⁿ⁺¹[nx+1,0,:] .= h[nx,0,:]
    # hⁿ⁺¹[0,ny+1,:] .= h[0,ny,:]
    # hⁿ⁺¹[nx+1,ny+1,:] .= h[nx,ny+1,:]

end

function cuda_wrap_kernel!(kernel::Function,
                           h, hⁿ⁺¹,
                           n1, n2, boundaryLength,
                           bc)

    ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
    iy = (blockIdx().y-1)*blockDim().y  + threadIdx().y
    if (0 < ix < n1+1 && 0 < iy < n2+1)
        kernel(ix, iy, h, hⁿ⁺¹, boundaryLength, bc)
    end

    return nothing
end


function boundaryloop!(kernel::Function,
                   h::CuArray,
                   hⁿ⁺¹::CuArray,
                   n1, n2, boundaryLength, bc)
    # Grid loop on GPU for CUarrays (CUDA)

    #ths = (8,8,4)
    ths = (8,8)
    nb1 = Int(ceil(n1/ths[1]))
    nb2 = Int(ceil(n2/ths[2]))
    bls = (nb1,nb2)
    @cuda threads=ths blocks=bls cuda_wrap_kernel!(kernel, h, hⁿ⁺¹, n1, n2, boundaryLength, bc)
end


function boundaryloop!(kernel::Function,
                   h::OffsetArray,
                   hⁿ⁺¹::OffsetArray,
                   n1, n2, boundaryLength, bc)
    # Grid loop for normal arrays on CPU

    for j = 1:n2
        for i = 1:n1
            kernel(i, j, h, hⁿ⁺¹, boundaryLength, bc)
        end
    end
end


function westBoundary!(j, k, h, hⁿ⁺¹, boundaryLength, bcval)
    h[0,j,k] = bcval
    hⁿ⁺¹[0,j,k] = bcval
end

function eastBoundary!(j, k, h, hⁿ⁺¹, boundaryLength, bcval)
    h[boundaryLength,j,k] = bcval
    hⁿ⁺¹[boundaryLength,j,k] = bcval
end

function southBoundary!(i, k, h, hⁿ⁺¹, boundaryLength, bcval)
    h[i,0,k] = h[i,1,k]
    hⁿ⁺¹[i,0,k] = h[i,1,k]
end

function northBoundary!(i, k, h, hⁿ⁺¹, boundaryLength, bcval)
    h[i,boundaryLength,k] = h[i,boundaryLength-1,k]
    hⁿ⁺¹[i,boundaryLength,k] = h[i,boundaryLength-1,k]
end

function bottomBoundary!(i, j, h, hⁿ⁺¹, boundaryLength, bcval)
    h[i,j,0] = h[i,j,1]
    hⁿ⁺¹[i,j,0] = h[i,j,1]
end

function topBoundary!(i, j, h, hⁿ⁺¹, boundaryLength, bcval)
    h[i,j,boundaryLength] = h[i,j,boundaryLength-1]
    hⁿ⁺¹[i,j,boundaryLength] = h[i,j,boundaryLength-1]
end
