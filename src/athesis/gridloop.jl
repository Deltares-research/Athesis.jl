using CUDAdrv, CUDAnative
using CuArrays

# This is the present data storage:
# grid       = (nx, ny, nz, Δx, Δy, Δz, x, y, z)
# state      = (h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
# model      = externals
# parameters = (K, Δt, tend)
# externals  = source
# source     = (i_src, j_src, k_src, n_src, externals)


function cuda_wrap_kernel!(kernel::Function,
                           source,
                           h,
                           u,
                           v,
                           w,
                           hⁿ⁺¹,
                           uⁿ⁺¹,
                           vⁿ⁺¹,
                           wⁿ⁺¹,
                           Δx::AbstractFloat,
                           Δy::AbstractFloat,
                           Δz::AbstractFloat,
                           nx,
                           ny,
                           nz,
                           Δt::AbstractFloat,
                           K,
                           SS::AbstractFloat)

    # 3D implementation
    ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
    iy = (blockIdx().y-1)*blockDim().y  + threadIdx().y
    iz = (blockIdx().z-1)*blockDim().z  + threadIdx().z
    if (0 < ix < nx+1 && 0 < iy < ny+1 && 0 < iz < nz+1)
        kernel(ix, iy, iz,
               source,
               h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹,
               Δx, Δy, Δz, Δt, K, SS)
    end

    return nothing
end


function gridloop!(kernel::Function,
                   source::CuArray,
                   state::State,
                   grid::Grid,
                   parameters::Parameters,
                   time_data::Time_data)
    # Grid loop on GPU for CUarrays (CUDA)

    # Unpack required variables
    nx    = grid.nx
    ny    = grid.ny
    nz    = grid.nz
    Δx    = grid.Δx
    Δy    = grid.Δy
    Δz    = grid.Δz
    Δt    = time_data.Δt
    h     = state.h
    u     = state.u
    v     = state.v
    w     = state.w
    hⁿ⁺¹  = state.hⁿ⁺¹
    uⁿ⁺¹  = state.uⁿ⁺¹
    vⁿ⁺¹  = state.vⁿ⁺¹
    wⁿ⁺¹  = state.wⁿ⁺¹
    K     = parameters.K
    SS    = parameters.specific_storage

    ths = (8,8,4)
    nbx = Int(ceil(grid.nx/ths[1]))
    nby = Int(ceil(grid.ny/ths[2]))
    nbz = Int(ceil(grid.nz/ths[3]))
    bls = (nbx,nby,nbz)
    @cuda threads=ths blocks=bls cuda_wrap_kernel!(
                            kernel, source,
                            h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹,
                            Δx, Δy, Δz,
                            nx, ny, nz, Δt,
                            K, SS)
end


function gridloop!(kernel::Function,
                   source::Array,
                   state::State,
                   grid::Grid,
                   parameters::Parameters,
                   time_data::Time_data)
    # Grid loop for normal arrays on CPU

    # Unpack required variables
    nx    = grid.nx
    ny    = grid.ny
    nz    = grid.nz
    Δx    = grid.Δx
    Δy    = grid.Δy
    Δz    = grid.Δz
    Δt    = time_data.Δt
    h     = state.h
    u     = state.u
    v     = state.v
    w     = state.w
    hⁿ⁺¹  = state.hⁿ⁺¹
    uⁿ⁺¹  = state.uⁿ⁺¹
    vⁿ⁺¹  = state.vⁿ⁺¹
    wⁿ⁺¹  = state.wⁿ⁺¹
    K     = parameters.K
    SS    = parameters.specific_storage

    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                kernel(i, j, k,
                       source,
                       h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹,
                       Δx, Δy, Δz, Δt, K, SS)
            end
        end
    end
end
