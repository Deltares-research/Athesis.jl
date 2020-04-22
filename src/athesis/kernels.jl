# kernels.jl

include("operators.jl")

include("kernels_pressure_equation.jl")
include("kernels_darcy_equation.jl")



# @inline function advection((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
#     # u*du/dx
#     uxm = (h[i,j]-h[i-1,j])/dx
#     uxp = (h[i+1,j]-h[i,j])/dx
#     am  = 0.5*min(u[i,j]+u[i+1,j],0.0)
#     ap  = 0.5*max(u[i-1,j]+u[i,j],0.0)
#     advx = ap*uxm + am*uxp
#
#     # v*du/dy
#     uym = (h[i,j]-h[i,j-1])/dy
#     uyp = (h[i,j+1]-h[i,j])/dy
#     am  = 0.5*min(v[i,j]+v[i,j+1],0.0)
#     ap  = 0.5*max(v[i,j-1]+v[i,j],0.0)
#     advy = ap*uym + am*uyp
#
#     #adv = 0.5*( (u[i,j]+u[i-1,j])*(u[i,j]-u[i-1,j])/dx + (v[i,j]+v[i,j-1])*(u[i,j]-u[i,j-1])/dy)
#     @inbounds result[i,j] -= dt*(advx+advy)
# end
#
# @inline function diffusion((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
#     @inbounds result[i,j] += dt*K*( (h[i+1,j] + h[i-1,j] - 2.0*h[i,j])/(dx*dx) + (h[i,j+1]+h[i,j-1] - 2.0*h[i,j])/(dy*dy) )
# end
#
# @inline function kernel_fused((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
#     result[i,j] = h[i,j]
#     # u*du/dx
#     uxm = (h[i,j]-h[i-1,j])/dx
#     uxp = (h[i+1,j]-h[i,j])/dx
#     am  = 0.5*min(u[i,j]+u[i+1,j],0.0)
#     ap  = 0.5*max(u[i-1,j]+u[i,j],0.0)
#     advx = ap*uxm + am*uxp
#
#     # v*du/dy
#     uym = (h[i,j]-h[i,j-1])/dy
#     uyp = (h[i,j+1]-h[i,j])/dy
#     am  = 0.5*min(v[i,j]+v[i,j+1],0.0)
#     ap  = 0.5*max(v[i,j-1]+v[i,j],0.0)
#     advy = ap*uym + am*uyp
#
#     #adv = 0.5*( (u[i,j]+u[i-1,j])*(u[i,j]-u[i-1,j])/dx + (v[i,j]+v[i,j-1])*(u[i,j]-u[i,j-1])/dy)
#     @inbounds result[i,j] -= dt*(advx+advy)
#     @inbounds result[i,j] += dt*K*( (h[i+1,j] + h[i-1,j] - 2.0*h[i,j])/(dx*dx) + (h[i,j+1]+h[i,j-1] - 2.0*h[i,j])/(dy*dy) )
# end
#
# @inline function kernel_composed((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
#     result[i,j] = h[i,j]
#     diffusion((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
#     advection((i,j), result, (h,u,v), (dx,dy,dt,K,nx,ny))
# end
