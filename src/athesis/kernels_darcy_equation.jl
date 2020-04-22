# kernels_darcy_equation.jl

include("operators.jl")


############## CPU kernels

@inline function darcy_kernel(i, j, source::Array, h::Array, u::Array, v::Array, hⁿ⁺¹::Array, uⁿ⁺¹::Array, vⁿ⁺¹::Array, Δx::Float64, Δy::Float64, Δt::Float64, K::Array)
    # 2D implementation (2DH or 2DV)

    uⁿ⁺¹[i,j] = -K[i,j]*∇2dᶜ(h,i,j,(Δx,Δy))[1]    # u
    vⁿ⁺¹[i,j] = -K[i,j]*∇2dᶜ(h,i,j,(Δx,Δy))[2]    # v or w
end

# @inline function darcy_kernel(i, k, source::Array{Float64,2}, h::Array{Float64,2}, u::Array{Float64,2}, w::Array{Float64,2}, hⁿ⁺¹::Array{Float64,2}, uⁿ⁺¹::Array{Float64,2}, wⁿ⁺¹::Array{Float64,2}, Δx::Float64, Δz::Float64, Δt::Float64, K::Array{Float64,2})
#     # 2DV implementation
#
#     uⁿ⁺¹[i,k] = -K[i,k]*∇2dᶜ(h,i,k,(Δx,Δz))[1]    # u
#     wⁿ⁺¹[i,k] = -K[i,k]*∇2dᶜ(h,i,k,(Δx,Δz))[2]    # w
# end

@inline function darcy_kernel(i, j, k, source::Array{Float64,3}, h::Array{Float64,3}, u::Array{Float64,3}, v::Array{Float64,3}, w::Array{Float64,3}, hⁿ⁺¹::Array{Float64,3}, uⁿ⁺¹::Array{Float64,3}, vⁿ⁺¹::Array{Float64,3}, wⁿ⁺¹::Array{Float64,3}, Δx::Float64, Δy::Float64, Δz::Float64, Δt::Float64, K::Array{Float64,3})
    # 3D implementation

    uⁿ⁺¹[i,j,k] = -K[i,j,k]*∇3dᶜ(h,i,j,k,(Δx,Δy,Δz))[1]    # u
    vⁿ⁺¹[i,j,k] = -K[i,j,k]*∇3dᶜ(h,i,j,k,(Δx,Δy,Δz))[2]    # v
    wⁿ⁺¹[i,j,k] = -K[i,j,k]*∇3dᶜ(h,i,j,k,(Δx,Δy,Δz))[3]    # w
end

############### GPU implementations

@inline function darcy_kernel(i, j, source::CuArray, h::CuArray, u::CuArray, v::CuArray, hⁿ⁺¹::CuArray, uⁿ⁺¹::CuArray, vⁿ⁺¹::CuArray, Δx::Float64, Δy::Float64, Δt::Float64, K::CuArray)
    # 2D implementation (2DH or 2DV)

    uⁿ⁺¹[i,j] = -K[i,j]*∇2dᶜ(h,i,j,(Δx,Δy))[1]    # u
    vⁿ⁺¹[i,j] = -K[i,j]*∇2dᶜ(h,i,j,(Δx,Δy))[2]    # v or w
end

# @inline function darcy_kernel(i, k, source::Array{Float64,2}, h::Array{Float64,2}, u::Array{Float64,2}, w::Array{Float64,2}, hⁿ⁺¹::Array{Float64,2}, uⁿ⁺¹::Array{Float64,2}, wⁿ⁺¹::Array{Float64,2}, Δx::Float64, Δz::Float64, Δt::Float64, K::Array{Float64,2})
#     # 2DV implementation
#
#     uⁿ⁺¹[i,k] = -K[i,k]*∇2dᶜ(h,i,k,(Δx,Δz))[1]    # u
#     wⁿ⁺¹[i,k] = -K[i,k]*∇2dᶜ(h,i,k,(Δx,Δz))[2]    # w
# end

@inline function darcy_kernel(i, j, k, source::CuArray{Float64,3}, h::CuArray{Float64,3}, u::CuArray{Float64,3}, v::CuArray{Float64,3}, w::CuArray{Float64,3}, hⁿ⁺¹::CuArray{Float64,3}, uⁿ⁺¹::CuArray{Float64,3}, vⁿ⁺¹::CuArray{Float64,3}, wⁿ⁺¹::CuArray{Float64,3}, Δx::Float64, Δy::Float64, Δz::Float64, Δt::Float64, K::CuArray{Float64,3})
    # 3D implementation

    uⁿ⁺¹[i,j,k] = -K[i,j,k]*∇3dᶜ(h,i,j,k,(Δx,Δy,Δz))[1]    # u
    vⁿ⁺¹[i,j,k] = -K[i,j,k]*∇3dᶜ(h,i,j,k,(Δx,Δy,Δz))[2]    # v
    wⁿ⁺¹[i,j,k] = -K[i,j,k]*∇3dᶜ(h,i,j,k,(Δx,Δy,Δz))[3]    # w
end
