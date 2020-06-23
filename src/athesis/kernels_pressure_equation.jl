# kernels_pressure_equation.jl

include("operators.jl")

########## CPU kernels
@inline function pressure_kernel!(i, j, source, h, u, v, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, Δx::Float64, Δy::Float64, Δt::Float64, K)
    # 2D implementation (2DH or 2DV)
    # for now cell centered
    F = K[i,j] * ( (h[i+1,j] + h[i-1,j] - 2.0*h[i,j])/(Δx*Δx) + (h[i,j+1]+h[i,j-1] - 2.0*h[i,j])/(Δy*Δy) ) + source[i,j]/(Δx*Δy)
    #F         = -K[i,j]*div2dᶜ(h,i,j,(Δx,Δy)) + source[i,j]
    hⁿ⁺¹[i,j] = h[i,j] + Δt*F
end

# @inline function pressure_kernel(i, k, source::Array{Float64,2}, h::Array{Float64,2}, u::Array{Float64,2}, w::Array{Float64,2}, hⁿ⁺¹::Array{Float64,2}, uⁿ⁺¹::Array{Float64,2}, wⁿ⁺¹::Array{Float64,2}, Δx::Float64, Δz::Float64, Δt::Float64, K::Array{Float64,2})
#     # 2DV implementation
#     # for now cell centered
#     F = -K[i,j] * ( (h[i+1,j] + h[i-1,j] - 2.0*h[i,j])/(Δx*Δx) + (h[i,j+1]+h[i,j-1] - 2.0*h[i,j])/(Δy*Δy) ) + source[i,k]
#     #F         = -K[i,k]*div2dᶜ(h,i,k,(Δx,Δz)) + source[i,k]
#     hⁿ⁺¹[i,k] = h[i,k] + Δt*F
# end

@inline function pressure_kernel!(i, j, k, source, h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹, Δx::Float64, Δy::Float64, Δz::Float64, Δt::Float64, K)
    # 3D implementation
    # for now cell centered
    F = K[i,j,k] * (
        (h[i+1,j,k] + h[i-1,j,k] - 2.0*h[i,j,k])/(Δx*Δx) +
        (h[i,j+1,k] + h[i,j-1,k] - 2.0*h[i,j,k])/(Δy*Δy) +
        (h[i,j,k+1] + h[i,j,k-1] - 2.0*h[i,j,k])/(Δz*Δz)
        ) + source[i,j,k]/(Δx*Δy)
    #F  = -K[i,j,k]*div3dᶜ(h,i,j,k, (Δx,Δy,Δz)) + source[i,j,k]
    hⁿ⁺¹[i,j,k] = h[i,j,k] + Δt*F
end

############# GPU kernels

# @inline function pressure_kernel(i, j, source::CuArray, h::CuArray, u::CuArray, v::CuArray, hⁿ⁺¹::CuArray, uⁿ⁺¹::CuArray, vⁿ⁺¹::CuArray, Δx::Float64, Δy::Float64, Δt::Float64, K::CuArray)
#     # 2D implementation (2DH or 2DV)
#     # for now cell centered
#     F = K[i,j] * ( (h[i+1,j] + h[i-1,j] - 2.0*h[i,j])/(Δx*Δx) + (h[i,j+1]+h[i,j-1] - 2.0*h[i,j])/(Δy*Δy) ) + source[i,j]/(Δx*Δy)
#     #F         = -K[i,j]*div2dᶜ(h,i,j,(Δx,Δy)) + source[i,j]
#     SS = 0.5
#     hⁿ⁺¹[i,j] = h[i,j] + Δt*F/SS
# end
#
# # @inline function pressure_kernel(i, k, source::CuArray{Float64,2}, h::CuArray{Float64,2}, u::CuArray{Float64,2}, w::CuArray{Float64,2}, hⁿ⁺¹::CuArray{Float64,2}, uⁿ⁺¹::CuArray{Float64,2}, wⁿ⁺¹::CuArray{Float64,2}, Δx::Float64, Δz::Float64, Δt::Float64, K::CuArray{Float64,2})
# #     # 2DV implementation
# #     # for now cell centered
# #     F = -K[i,j] * ( (h[i+1,j] + h[i-1,j] - 2.0*h[i,j])/(Δx*Δx) + (h[i,j+1]+h[i,j-1] - 2.0*h[i,j])/(Δy*Δy) ) + source[i,k]
# #     #F         = -K[i,k]*div2dᶜ(h,i,k,(Δx,Δz)) + source[i,k]
# #     hⁿ⁺¹[i,k] = h[i,k] + Δt*F
# # end
#
# @inline function pressure_kernel(i, j, k, source::CuArray{Float64,3}, h::CuArray{Float64,3}, u::CuArray{Float64,3}, v::CuArray{Float64,3}, w::CuArray{Float64,3}, hⁿ⁺¹::CuArray{Float64,3}, uⁿ⁺¹::CuArray{Float64,3}, vⁿ⁺¹::CuArray{Float64,3}, wⁿ⁺¹::CuArray{Float64,3}, Δx::Float64, Δy::Float64, Δz::Float64, Δt::Float64, K::CuArray{Float64,3})
#     # 3D implementation
#     # for now cell centered
#     F = -K[i,j,k] * (
#         (h[i+1,j,k] + h[i-1,j,k] - 2.0*h[i,j,k])/(Δx*Δx) +
#         (h[i,j+1,k] + h[i,j-1,k] - 2.0*h[i,j,k])/(Δy*Δy) +
#         (h[i,j,k+1] + h[i,j,k-1] - 2.0*h[i,j,k])/(Δz*Δz)
#         ) + source[i,j,k]
#     #F  = -K[i,j,k]*div3dᶜ(h,i,j,k, (Δx,Δy,Δz)) + source[i,j,k]
#     hⁿ⁺¹[i,j,k] = h[i,j,k] + Δt*F
# end
