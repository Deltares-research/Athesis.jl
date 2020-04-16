using CUDAdrv, CUDAnative
using CuArrays
using TimerOutputs

function backendname(useCUDA)
    return useCUDA ? "GPU" : "CPU"
end

function run()

    to = TimerOutput()

    nx = 1000
    ny = 1000
    N = 1000000

    # Model parameters
    dx = 10.0
    dy = 10.0
    dt = 10.0
    K  = 1.0e-4

    println("Hit c for cuda...")
    c = readline()
    useCUDA = (c == "c")

    h = fill(5.0, N)
    # for i = 1:N
    #     h[i] = 10.0*i/N
    # end
    u = fill(2.0, N)
    v = fill(1.0, N)

    if useCUDA
        h = CuArray(h)
        u = CuArray(u)
        v = CuArray(v)
    end
    result = similar(h)
    state  = (h,u,v)
    params = (dx, dy, dt, K, nx, ny)

    # first separate kernels:
    fill!(result, 0.0)
    println("separate, in: ", h[1:10])
    @timeit to "separate "*backendname(useCUDA) begin
    gridloop(advection, result, state, params)
    gridloop(diffusion, result, state, params)
    end
    println("result:  ", result[1:10])

    # then merged into one:
    fill!(result, 0.0)
    println("fused, in: ", h[1:10])
    @timeit to "fused "*backendname(useCUDA) begin
    gridloop(kernel_fused, result, state, params)
    end
    println("result:  ", result[1:10])

    # then composed:
    fill!(result, 0.0)
    println("composed, in: ", h[1:10])
    @timeit to "composed "*backendname(useCUDA) begin
    gridloop(kernel_composed, result, state, params)
    end
    println("result:  ", result[1:10])

    print_timer(to)
end

@inline function advection((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
    # u*du/dx
    uxm = (h[i,j]-h[i-1,j])/dx
    uxp = (h[i+1,j]-h[i,j])/dx
    am  = 0.5*min(u[i,j]+u[i+1,j],0.0)
    ap  = 0.5*max(u[i-1,j]+u[i,j],0.0)
    advx = ap*uxm + am*uxp

    # v*du/dy
    uym = (h[i,j]-h[i,j-1])/dy
    uyp = (h[i,j+1]-h[i,j])/dy
    am  = 0.5*min(v[i,j]+v[i,j+1],0.0)
    ap  = 0.5*max(v[i,j-1]+v[i,j],0.0)
    advy = ap*uym + am*uyp

    #adv = 0.5*( (u[i,j]+u[i-1,j])*(u[i,j]-u[i-1,j])/dx + (v[i,j]+v[i,j-1])*(u[i,j]-u[i,j-1])/dy)
    @inbounds result[i,j] -= dt*(advx+advy)
end

@inline function diffusion((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
    @inbounds result[i,j] += dt*K*( (h[i+1,j] + h[i-1,j] - 2.0*h[i,j])/(dx*dx) + (h[i,j+1]+h[i,j-1] - 2.0*h[i,j])/(dy*dy) )
end

@inline function kernel_fused((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
    result[i,j] = h[i,j]
    # u*du/dx
    uxm = (h[i,j]-h[i-1,j])/dx
    uxp = (h[i+1,j]-h[i,j])/dx
    am  = 0.5*min(u[i,j]+u[i+1,j],0.0)
    ap  = 0.5*max(u[i-1,j]+u[i,j],0.0)
    advx = ap*uxm + am*uxp

    # v*du/dy
    uym = (h[i,j]-h[i,j-1])/dy
    uyp = (h[i,j+1]-h[i,j])/dy
    am  = 0.5*min(v[i,j]+v[i,j+1],0.0)
    ap  = 0.5*max(v[i,j-1]+v[i,j],0.0)
    advy = ap*uym + am*uyp

    #adv = 0.5*( (u[i,j]+u[i-1,j])*(u[i,j]-u[i-1,j])/dx + (v[i,j]+v[i,j-1])*(u[i,j]-u[i,j-1])/dy)
    @inbounds result[i,j] -= dt*(advx+advy)
    @inbounds result[i,j] += dt*K*( (h[i+1,j] + h[i-1,j] - 2.0*h[i,j])/(dx*dx) + (h[i,j+1]+h[i,j-1] - 2.0*h[i,j])/(dy*dy) )
end

@inline function kernel_composed((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
    result[i,j] = h[i,j]
    diffusion((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
    advection((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
end

function cuda_wrap_kernel(kernel::Function, result::AbstractArray{Float64}, state, args)
    nx = args[5]
    ny = args[6]

    ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
    iy = (blockIdx().y-1)*blockDim().y  + threadIdx().y
    if (ix <= nx && iy && ny)
        kernel((ix,iy), result, state, args)
    end

    iy = blockDim().y * (blockIdx().y - 1) + threadIdx().y;
    idx = size(u,1)*(iy-1) + ix
    if ix < size(u,1) + 1 && iy < size(u,2) + 1
        @inbounds w[idx] = abs(v[ix,iy] - u[ix,iy])
    end

    return nothing
end

# function gridloop(kernel::Function, result::CuArray{Float64}, state, args)
#     # Grid loop for normal arrays on CPU for CUarrays (CUDA)
#     nx = args[5]
#     ny = args[6]
#     result = reshape(result, (nx,ny))
#     h = state[1]
#     u = state[2]
#     v = state[3]
#     h = reshape(h, (nx,ny))
#     u = reshape(u, (nx,ny))
#     v = reshape(v, (nx,ny))
#     state = (h,u,v)
#
#     ths = 256
#     bls = Int(ceil(length(result) / ths))
#     @cuda threads=ths blocks=bls cuda_wrap_kernel(kernel, result, state, args)
# end
#
# function gridloop(kernel::Function, result::Array{Float64}, state, args)
#     # Grid loop for normal arrays on CPU for normal arrays
#     nx = args[5]
#     ny = args[6]
#     result = reshape(result, (nx,ny))
#     h = state[1]
#     u = state[2]
#     v = state[3]
#     h = reshape(h, (nx,ny))
#     u = reshape(u, (nx,ny))
#     v = reshape(v, (nx,ny))
#     state = (h,u,v)
#
#     for j = 2:ny-1
#         for i = 2:nx-1
#             kernel((i,j), result, state, args)
#         end
#     end
#     result = vec(result)
#     h = vec(h)
#     u = vec(u)
#     v = vec(v)
# end

function gridloop(kernel::Function, result::AbstractArray{Float64}, args)
    nx = args[5]
    ny = args[6]
    result = reshape(result, (nx,ny))
    h = state[1]
    u = state[2]
    v = state[3]
    h = reshape(h, (nx,ny))
    u = reshape(u, (nx,ny))
    v = reshape(v, (nx,ny))
    state = (h,u,v)

    if result isa CuArray
        ths = 256
        bls = Int(ceil(length(result) / ths))
        @cuda threads=ths blocks=bls cuda_wrap_kernel(kernel, result, args)
    else
        for i = 1 : length(result)
            kernel(i, result, args)
        end
    end
    result = vec(result)
    h = vec(h)
    u = vec(u)
    v = vec(v)
end
