# grids.jl

mutable struct Grid{T}
    #grid_type::String
    nx::Int64
    ny::Int64
    nz::Int64
    Δx::Float64
    Δy::Float64
    Δz::Float64
    x::T
    y::T
    z::T
end

function grid_coords(nx, ny, nz, Δx, Δy, Δz, useCUDA)
    # Create 3D structured grid
    # For now cell centered
    x = Array{Float64,1}(undef,nx+2)
    for n = 1:nx+2
        x[n] = (n-1.5)*Δx
    end

    y = Array{Float64,1}(undef,ny+2)
    for n = 1:ny+2
        y[n] = (n-1.5)*Δy
    end

    z = Array{Float64,1}(undef,nz+2)
    for n = 1:nz+2
        z[n] = (n-1.5)*Δz
    end

    if useCUDA
        # Convert to CUDA Arrays
        x = CuArray(x)
        y = CuArray(y)
        z = CuArray(z)
    end
    return x, y, z
end


# function sigma_grid!(grid, state)
#     nx  = grid[1]
#     nz  = grid[2]
#     Δx  = grid.Δx
#     Δz  = grid.Δz
#     xc  = grid.xc
#     zc  = grid.zc
#     xu  = grid.xu
#     zu  = grid.zu
#     xw  = grid.xw
#     zw  = grid.zw
#
#     dzs = state.Δzₛ
#     dzu = state.Δzᵤ
#
#     # Co-ordinates of the cell centres
#     zc[:,1] = 0.5*dzs[:,1]
#     for i = 1:nx
#         #for j = 1:ny
#         for k = 2:nz
#             xc[i,k] = (i-0.5)*Δx
#             zc[i,k] = zc[i,k-1] + 0.5*(dzs[i,k-1]+dzs[i,k])
#         end
#     end
#
#     # Co-ordinates of the horizontal velocity points
#     zu[:,1] = 0.5*dzu[:,1]
#     for i = 1:nx+1
#         #for j = 1:ny
#         for k = 2:nz
#             xu[i,k] = (i-1)*Δx
#             zu[i,k] = zu[i,k-1] + 0.5*(dzu[i,k-1]+dzu[i,k])
#         end
#     end
#
#     # Co-ordinates of the vertical velocity points
#     zw[:,1] .= 0.0
#     for i = 1:nx
#         #for j = 1:ny
#         for k = 2:nz+1
#             xw[i,k] = (i-0.5)*Δx
#             zw[i,k] = zw[i,k-1] + dzs[i,k-1]
#         end
#     end
#
# end
