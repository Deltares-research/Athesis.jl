#copy_state.jl

using CuArrays

function new2old!(state::State,bc::BoundaryConditions)
    h    = state.h
    u    = state.u
    v    = state.v
    w    = state.w
    hⁿ⁺¹ = state.hⁿ⁺¹
    uⁿ⁺¹ = state.uⁿ⁺¹
    vⁿ⁺¹ = state.vⁿ⁺¹
    wⁿ⁺¹ = state.wⁿ⁺¹

    bcpres = bc.bc_pressure

    copy_new2old!(bcpres, h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
end

function copy_new2old!(bcpres::Array, h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)

    #println("CPU copy")
    for k = 0:size(h)[3]-1
        for j = 0:size(h)[2]-1
            for i = 0:size(h)[1]-1
                copy_kernel!(i, j, k,
                       h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
            end
        end
    end
end

function copy_new2old!(bcpres::CuArray, h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)

    #println("GPU copy")
    ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
    iy = (blockIdx().y-1)*blockDim().y  + threadIdx().y
    iz = (blockIdx().z-1)*blockDim().z  + threadIdx().z

    ths = (8,8,4)
    nbx = Int(ceil(size(h)[1]/ths[1]))
    nby = Int(ceil(size(h)[2]/ths[2]))
    nbz = Int(ceil(size(h)[3]/ths[3]))
    bls = (nbx,nby,nbz)

    @cuda threads=ths blocks=bls cuda_wrap_copy!(
                            h, u, v, w,
                            hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)

end

function cuda_wrap_copy!(h, u, v, w,
                         hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)

    ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
    iy = (blockIdx().y-1)*blockDim().y  + threadIdx().y
    iz = (blockIdx().z-1)*blockDim().z  + threadIdx().z
    if (0 <= ix <= nx+1 && 0 <= iy <= ny+1 && 0 <= iz <= nz+1)
        copy_kernel!(ix, iy, iz,
               h, u, v, w,
               hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
    end

    return nothing
end

function copy_kernel!(i,j,k,
                      h, u, v, w,
                      hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹ )

    h[i,j,k] = hⁿ⁺¹[i,j,k]
    u[i,j,k] = uⁿ⁺¹[i,j,k]
    v[i,j,k] = vⁿ⁺¹[i,j,k]
    w[i,j,k] = wⁿ⁺¹[i,j,k]
end
