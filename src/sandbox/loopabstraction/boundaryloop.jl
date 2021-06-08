# Treatment of boundary conditions
# Using abstract loops on either CPU or GPU

boundaryloop!(westBoundary!, state.h, state.hⁿ⁺¹, ny, nz, nx + 1, bc.bc_pressure)

boundaryloop!(eastBoundary!, state.h, state.hⁿ⁺¹, ny, nz, nx + 1, bc.bc_pressure)

boundaryloop!(southBoundary!, state.h, state.hⁿ⁺¹, nx, nz, ny + 1, bc.bc_pressure)

boundaryloop!(northBoundary!, state.h, state.hⁿ⁺¹, nx, nz, ny + 1, bc.bc_pressure)

boundaryloop!(bottomBoundary!, state.h, state.hⁿ⁺¹, nx, ny, nz + 1, bc.bc_pressure)

boundaryloop!(topBoundary!, state.h, state.hⁿ⁺¹, nx, ny, nz + 1, bc.bc_pressure)

function cuda_wrap_kernel!(kernel::Function,
                           h, hⁿ⁺¹,
                           n1, n2, boundaryLength,
                           bc::CuArray)

    ix = (blockIdx().x - 1) * blockDim().x  + threadIdx().x
    iy = (blockIdx().y - 1) * blockDim().y  + threadIdx().y
    if (0 < ix < n1 + 1 && 0 < iy < n2 + 1)
        kernel(ix, iy, h, hⁿ⁺¹, boundaryLength, bc)
    end

    return nothing
end


function boundaryloop!(kernel::Function,
                   h,
                   hⁿ⁺¹,
                   n1, n2, boundaryLength,
                   bc::CuArray)
    # Grid loop on GPU for CUarrays (CUDA)

    # ths = (8,8,4)
    ths = (8, 8)
    nb1 = Int(ceil(n1 / ths[1]))
    nb2 = Int(ceil(n2 / ths[2]))
    bls = (nb1, nb2)
    @cuda threads = ths blocks = bls cuda_wrap_kernel!(kernel, h, hⁿ⁺¹, n1, n2, boundaryLength, bc)
end


function boundaryloop!(kernel::Function,
                   h,
                   hⁿ⁺¹,
                   n1, n2, boundaryLength,
                   bc::Array)
    # Grid loop for normal arrays on CPU

    for j = 1:n2
        for i = 1:n1
            kernel(i, j, h, hⁿ⁺¹, boundaryLength, bc)
        end
    end
end


function westBoundary!(j, k, h, hⁿ⁺¹, boundaryLength, bc)
    h[0,j,k] = bc[1]
    hⁿ⁺¹[0,j,k] = bc[1]
end

function eastBoundary!(j, k, h, hⁿ⁺¹, boundaryLength, bc)
    h[boundaryLength,j,k] = bc[2]
    hⁿ⁺¹[boundaryLength,j,k] = bc[2]
end

function southBoundary!(i, k, h, hⁿ⁺¹, boundaryLength, bc)
    h[i,0,k] = h[i,1,k]
    hⁿ⁺¹[i,0,k] = h[i,1,k]
end

function northBoundary!(i, k, h, hⁿ⁺¹, boundaryLength, bc)
    h[i,boundaryLength,k] = h[i,boundaryLength - 1,k]
    hⁿ⁺¹[i,boundaryLength,k] = h[i,boundaryLength - 1,k]
end

function bottomBoundary!(i, j, h, hⁿ⁺¹, boundaryLength, bc)
    h[i,j,0] = h[i,j,1]
    hⁿ⁺¹[i,j,0] = h[i,j,1]
end

function topBoundary!(i, j, h, hⁿ⁺¹, boundaryLength, bc)
    h[i,j,boundaryLength] = h[i,j,boundaryLength - 1]
    hⁿ⁺¹[i,j,boundaryLength] = h[i,j,boundaryLength - 1]
end
