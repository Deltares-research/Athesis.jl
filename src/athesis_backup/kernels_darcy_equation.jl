# kernels_darcy_equation.jl

include("operators.jl")

############## Darcy kernel (CPU and GPU compatible)

@inline function darcy_kernel!(i, j, k, source, h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹, Δx::Float64, Δy::Float64, Δz::Float64, Δt::Float64, K, SS)
    uⁿ⁺¹[i,j,k] = -K[i,j,k]*∇3dᶜ(h,i,j,k,(Δx,Δy,Δz))[1]    # u
    vⁿ⁺¹[i,j,k] = -K[i,j,k]*∇3dᶜ(h,i,j,k,(Δx,Δy,Δz))[2]    # v
    wⁿ⁺¹[i,j,k] = -K[i,j,k]*∇3dᶜ(h,i,j,k,(Δx,Δy,Δz))[3]    # w
end
