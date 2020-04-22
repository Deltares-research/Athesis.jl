# Generic operators (e.g. advection and diffusion, explicit/implicit)
# using building blocks:
# - difference operator
# - gradient
# - divergence
# - laplacian

# Central and upwind types
#struct Central  end; central  = Central()
#struct Upwind   end; upwind   = Upwind()
#struct Downwind end; downwind = Downwind()

############### 2d difference operators ##################
# Difference operators in x- and y-direction:
# Central-x
@inline δfδx_2dᶜ(f,i,j,Δx) = (f[i+1,j] - f[i-1,j])/Δx
# Upwind-x
@inline δfδx_2dᵘ(f,i,j,Δx) = (f[i,j] - f[i-1,j])/Δx
# Downwind-x
@inline δfδx_2dᵈ(f,i,j,Δx) = (f[i+1,j] - f[i,j])/Δx
# Central-y
@inline δfδy_2dᶜ(f,i,j,Δy) = (f[i,j+1] - f[i,j-1])/Δy
# Upwind-y
@inline δfδy_2dᵘ(f,i,j,Δy) = (f[i,j] - f[i,j-1])/Δy
# Downwind-y
@inline δfδy_2dᵈ(f,i,j,Δy) = (f[i,j+1] - f[i,j])/Δy

############### 3d difference operators ##################
# Difference operators in x- and y-direction:
# Central-x
@inline δfδx_3dᶜ(f,i,j,k,Δx)  = (f[i+1,j,k] - f[i-1,j,k])/Δx
# Upwind-x
@inline δfδx_3dᵘ(f,i,j,k,Δx)   = (f[i,j,k] - f[i-1,j,k])/Δx
# Downwind-x
@inline δfδx_3dᵈ(f,i,j,k,Δx) = (f[i+1,j,k] - f[i,j,k])/Δx
# Central-y
@inline δfδy_3dᶜ(f,i,j,k,Δy)  = (f[i,j+1,k] - f[i,j-1,k])/Δy
# Upwind-y
@inline δfδy_3dᵘ(f,i,j,k,Δy)   = (f[i,j,k] - f[i,j-1,k])/Δy
# Downwind-y
@inline δfδy_3dᵈ(f,i,j,k,Δy) = (f[i,j+1,k] - f[i,j,k])/Δy
# Central-z
@inline δfδz_3dᶜ(f,i,j,k,Δz)  = (f[i,j,k+1] - f[i,j,k-1])/Δz
# Upwind-z
@inline δfδz_3dᵘ(f,i,j,k,Δz)   = (f[i,j,k] - f[i,j,k-1])/Δz
# Downwind-z
@inline δfδz_3dᵈ(f,i,j,k,Δz) = (f[i,j,k+1] - f[i,j,k])/Δz

############### 2d gradient operators ##################
# central gradient
@inline ∇2dᶜ(f,i,j,(Δx,Δy)) = (δfδx_2dᶜ(f,i,j,Δx), δfδy_2dᶜ(f,i,j,Δy))
# upwind gradient
@inline ∇2dᵘ(f,i,j,(Δx,Δy)) = (δfδx_2dᵘ(f,i,j,Δx), δfδy_2dᵘ(f,i,j,Δy))
# downwind gradient
@inline ∇2dᵈ(f,i,j,(Δx,Δy)) = (δfδx_2dᵈ(f,i,j,Δx), δfδy_2dᵈ(f,i,j,Δy))

############### 3d gradient operators ##################
# central gradient
@inline ∇3dᶜ(f,i,j,k,(Δx,Δy,Δz)) = (δfδx_3dᶜ(f,i,j,k,Δx), δfδy_3dᶜ(f,i,j,k,Δy), δfδz_3dᶜ(f,i,j,k,Δz))
# upwind gradient
@inline ∇3dᵘ(f,i,j,k,(Δx,Δy,Δz)) = (δfδx_3dᵘ(f,i,j,k,Δx), δfδy_3dᵘ(f,i,j,k,Δy), δfδz_3dᵘ(f,i,j,k,Δz))
# downwind gradient
@inline ∇3dᵈ(f,i,j,k,(Δx,Δy,Δz)) = (δfδx_3dᵈ(f,i,j,k,Δx), δfδy_3dᵈ(f,i,j,k,Δy), δfδz_3dᵈ(f,i,j,k,Δz))

############### 2d divergence operators ##################
# central divergence
@inline div2dᶜ(f,i,j,(Δx,Δy)) = δfδx_2dᶜ(f,i,j,Δx) + δfδy_2dᶜ(f,i,j,Δy)
# upwind divergence
@inline div2dᵘ(f,i,j,(Δx,Δy)) = δfδx_2dᵘ(f,i,j,Δx) + δfδy_2dᵘ(f,i,j,Δy)
# downwind divergence
@inline div2dᵈ(f,i,j,(Δx,Δy)) = δfδx_2dᵈ(f,i,j,Δx) + δfδy_2dᵈ(f,i,j,Δy)

############### 3d divergence operators ##################
# central divergence
@inline div3dᶜ(f,i,j,k,(Δx,Δy,Δz)) = (δfδx_3dᶜ(f,i,j,k,Δx) + δfδy_3dᶜ(f,i,j,k,Δy) + δfδz_3dᶜ(f,i,j,k,Δz))
# upwind divergence
@inline div3dᵘ(f,i,j,k,(Δx,Δy,Δz)) = (δfδx_3dᵘ(f,i,j,k,Δx) + δfδy_3dᵘ(f,i,j,k,Δy) + δfδz_3dᵘ(f,i,j,k,Δz))
# downwind divergence
@inline div3dᵈ(f,i,j,k,(Δx,Δy,Δz)) = (δfδx_3dᵈ(f,i,j,k,Δx) + δfδy_3dᵈ(f,i,j,k,Δy) + δfδz_3dᵈ(f,i,j,k,Δz))


# Laplacian = div(gradient(f))
# central Laplacian
#@inline Δˣˣ(f,i,j,central)  = div(∇(f,i,j,central ),i,j,central)
#@inline Δˣˣ(f,i,j,upwind)   = div(∇(f,i,j,upwind  ),i,j,upwind)
#@inline Δˣˣ(f,i,j,downwind) = div(∇(f,i,j,downwind),i,j,downwind)
