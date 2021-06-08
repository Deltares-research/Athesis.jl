using KernelAbstractions

@kernel function darcyKernel!(uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹, h, Δx, Δy, Δz, K)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        uⁿ⁺¹[i,j,k] = -K[i,j,k] * ∇3dᶜ(h, i, j, k, (Δx, Δy, Δz))[1]    # u
        vⁿ⁺¹[i,j,k] = -K[i,j,k] * ∇3dᶜ(h, i, j, k, (Δx, Δy, Δz))[2]    # v
        wⁿ⁺¹[i,j,k] = -K[i,j,k] * ∇3dᶜ(h, i, j, k, (Δx, Δy, Δz))[3]    # w
    end
end
