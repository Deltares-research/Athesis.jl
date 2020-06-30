using CUDAdrv, CUDAnative
using CuArrays
using TimerOutputs

include("kernels.jl")
include("model.jl")
include("grids.jl")

# function backendname(useCUDA)
#     return useCUDA ? "GPU" : "CPU"
# end

# This is the present data storage:
# grid       = (ndim, nx, ny, Δx, Δy, x, y) or (ndim, nx, ny, nz, Δx, Δy, Δz, x, y, z)
# state      = (h, u, v, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹) or (h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
# model      = externals
# parameters = (K, Δt, tend)
# externals  = source
# source     = (i_src, j_src, n_src, externals) or (i_src, j_src, k_src, n_src, externals)

# function cuda_wrap_kernel!(kernel::Function,
#                            source::CuArray{Float64,2},
#                            h::CuArray{Float64,2},
#                            u::CuArray{Float64,2},
#                            v::CuArray{Float64,2},
#                            hⁿ⁺¹::CuArray{Float64,2},
#                            uⁿ⁺¹::CuArray{Float64,2},
#                            vⁿ⁺¹::CuArray{Float64,2},
#                            Δx::Float64,
#                            Δy::Float64,
#                            Δt::Float64,
#                            K::CuArray{Float64,2})
#
#     # 2DH implementation
#
#     ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
#     iy = (blockIdx().y-1)*blockDim().y  + threadIdx().y
#     if (1 < ix < nx && 1 < iy < ny)
#         kernel(ix, iy, source, h, u, v, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, Δx, Δy, Δt, K)
#     end
#
#     return nothing
# end

# function cuda_wrap_kernel!(kernel::Function,
#                            source::CuArray,
#                            h::CuArray,
#                            u::CuArray,
#                            w::CuArray,
#                            hⁿ⁺¹::CuArray,
#                            uⁿ⁺¹::CuArray,
#                            wⁿ⁺¹::CuArray,
#                            Δx::Float64,
#                            Δz::Float64,
#                            Δt::Float64,
#                            K::CuArray)
#
#     # 2DV implementation
#
#     ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
#     iy = (blockIdx().y-1)*blockDim().y  + threadIdx().y
#     if (1 < ix < nx && 1 < iy < nz)
#         kernel(ix, iz, source, h, u, w, hⁿ⁺¹, uⁿ⁺¹, wⁿ⁺¹, Δx, Δz, Δt, K)
#     end
#
#     return nothing
# end

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
                           Δx::Float64,
                           Δy::Float64,
                           Δz::Float64,
                           nx,
                           ny,
                           nz,
                           Δt::Float64,
                           K)

    # 3D implementation
    ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
    iy = (blockIdx().y-1)*blockDim().y  + threadIdx().y
    iz = (blockIdx().z-1)*blockDim().z  + threadIdx().z
    if (1 < ix < nx && 1 < iy < ny && 1 < iz < nz)
        kernel(ix,iy,iz,
               source,
               h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹,
               Δx, Δy, Δz, Δt, K)
    end

    return nothing
end

function gridloop!(kernel::Function,
                   source::CuArray,
                   state::State2dh, # , u::CuArray{Float64}, v::CuArray{Float64}, hⁿ⁺¹::CuArray{Float64}, uⁿ⁺¹::CuArray{Float64}, vⁿ⁺¹::CuArray{Float64},
                   grid::Grid2dh, #ndim::Int64, nx::Int64, ny::Int64, Δx::Float64, Δy::Float64, x::CuArray{Float64}, y::CuArray{Float64},
                   parameters::Parameters, #
                   time_data::Time_data) #Δt::Float64, tend::Float64)
    # Grid loop on CPU for CUarrays (CUDA)
    # 2DH implementation

    # Unpack required variables
    nx    = grid.nx
    ny    = grid.ny
    Δx    = grid.Δx
    Δy    = grid.Δy
    Δt    = time_data.Δt
    h     = state.h
    u     = state.u
    v     = state.v
    hⁿ⁺¹  = state.hⁿ⁺¹
    uⁿ⁺¹  = state.uⁿ⁺¹
    vⁿ⁺¹  = state.vⁿ⁺¹
    K     = parameters.K

    ths = 256
    bls = Int(ceil(length(h) / ths))
    @cuda threads=ths blocks=bls cuda_wrap_kernel!(kernel, source, h, u, v, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, Δx, Δy, Δt, K)
end

function gridloop!(kernel::Function,
                   source::CuArray,
                   state::State2dv, # , u::CuArray{Float64}, v::CuArray{Float64}, hⁿ⁺¹::CuArray{Float64}, uⁿ⁺¹::CuArray{Float64}, vⁿ⁺¹::CuArray{Float64},
                   grid::Grid2dv, #ndim::Int64, nx::Int64, ny::Int64, Δx::Float64, Δy::Float64, x::CuArray{Float64}, y::CuArray{Float64},
                   parameters::Parameters, #
                   time_data::Time_data) #Δt::Float64, tend::Float64)
    # Grid loop on CPU for CUarrays (CUDA)
    # 2DV implementation

    # Unpack required variables
    nx    = grid.nx
    nz    = grid.nz
    Δx    = grid.Δx
    Δz    = grid.Δz
    Δt    = time_data.Δt
    h     = state.h
    u     = state.u
    w     = state.w
    hⁿ⁺¹  = state.hⁿ⁺¹
    uⁿ⁺¹  = state.uⁿ⁺¹
    wⁿ⁺¹  = state.wⁿ⁺¹
    K     = parameters.K

    ths = 256
    bls = Int(ceil(length(h) / ths))
    @cuda threads=ths blocks=bls cuda_wrap_kernel!(kernel, source, h, u, w, hⁿ⁺¹, uⁿ⁺¹, wⁿ⁺¹, Δx, Δz, Δt, K)
end

function gridloop!(kernel::Function,
                   source::CuArray,
                   state::State3d, # , u::CuArray{Float64}, v::CuArray{Float64}, hⁿ⁺¹::CuArray{Float64}, uⁿ⁺¹::CuArray{Float64}, vⁿ⁺¹::CuArray{Float64},
                   grid::Grid3d, #ndim::Int64, nx::Int64, ny::Int64, Δx::Float64, Δy::Float64, x::CuArray{Float64}, y::CuArray{Float64},
                   parameters::Parameters, #
                   time_data::Time_data)
    # Grid loop on CPU for CUarrays (CUDA)
    # 3D implementation

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

    ths = (8,8,4)
    nbx = Int(ceil(grid.nx/ths[1]))
    nby = Int(ceil(grid.ny/ths[2]))
    nbz = Int(ceil(grid.nz/ths[3]))
    bls = (nbx,nby,nbz)
    @cuda threads=ths blocks=bls cuda_wrap_kernel!(kernel, source, h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹, Δx, Δy, Δz, nx, ny, nz, Δt, K)
end

function gridloop!(kernel::Function,
                   source::Array{Float64,2},
                   state::State2dh, # , u::CuArray{Float64}, v::CuArray{Float64}, hⁿ⁺¹::CuArray{Float64}, uⁿ⁺¹::CuArray{Float64}, vⁿ⁺¹::CuArray{Float64},
                   grid::Grid2dh, #ndim::Int64, nx::Int64, ny::Int64, Δx::Float64, Δy::Float64, x::CuArray{Float64}, y::CuArray{Float64},
                   parameters::Parameters, #
                   time_data::Time_data)
    # Grid loop for normal arrays on CPU for normal arrays

    # 2DH implementation

    # Unpack required variables
    nx    = grid.nx
    ny    = grid.ny
    Δx    = grid.Δx
    Δy    = grid.Δy
    Δt    = time_data.Δt
    h     = state.h
    u     = state.u
    v     = state.v
    hⁿ⁺¹  = state.hⁿ⁺¹
    uⁿ⁺¹  = state.uⁿ⁺¹
    vⁿ⁺¹  = state.vⁿ⁺¹
    K     = parameters.K

    for j = 2:ny-1
        for i = 2:nx-1
            kernel(i, j, source, h, u, v, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, Δx, Δy, Δt, K)
        end
    end
end

function gridloop!(kernel::Function,
                  source::Array{Float64,2},
                  state::State2dv, # , u::CuArray{Float64}, v::CuArray{Float64}, hⁿ⁺¹::CuArray{Float64}, uⁿ⁺¹::CuArray{Float64}, vⁿ⁺¹::CuArray{Float64},
                  grid::Grid2dv, #ndim::Int64, nx::Int64, ny::Int64, Δx::Float64, Δy::Float64, x::CuArray{Float64}, y::CuArray{Float64},
                  parameters::Parameters, #
                  time_data::Time_data)
    # Grid loop for normal arrays on CPU for normal arrays

    # 2DV implementation

    # Unpack required variables
    nx    = grid.nx
    nz    = grid.nz
    Δx    = grid.Δx
    Δz    = grid.Δz
    Δt    = time_data.Δt
    h     = state.h
    u     = state.u
    w     = state.w
    hⁿ⁺¹  = state.hⁿ⁺¹
    uⁿ⁺¹  = state.uⁿ⁺¹
    wⁿ⁺¹  = state.wⁿ⁺¹
    K     = parameters.K

    for k = 2:nz-1
        for i = 2:nx-1
            kernel(i, k, source, h, u, w, hⁿ⁺¹, uⁿ⁺¹, wⁿ⁺¹, Δx, Δz, Δt, K)
        end
    end
end

function gridloop!(kernel::Function,
                   source::Array{Float64,3},
                   state::State3d, # , u::CuArray{Float64}, v::CuArray{Float64}, hⁿ⁺¹::CuArray{Float64}, uⁿ⁺¹::CuArray{Float64}, vⁿ⁺¹::CuArray{Float64},
                   grid::Grid3d, #ndim::Int64, nx::Int64, ny::Int64, Δx::Float64, Δy::Float64, x::CuArray{Float64}, y::CuArray{Float64},
                   parameters::Parameters, #
                   time_data::Time_data)
    # Grid loop for normal arrays on CPU for normal arrays

    # 3D implementation

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

    for k = 2:nz-1
        for j = 2:ny-1
            for i = 2:nx-1
                kernel(i,j,k,
                       source,
                       h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹,
                       Δx, Δy, Δz, Δt, K)
            end
        end
    end
end
