# grids.jl

using Adapt

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

function grid_coords(n1, n2, n3, Δ1, Δ2, Δ3, useCUDA)
    # Create 3D structured grid
    # For now cell centered
    coord1 = Array{Float64,1}(undef,n1)
    for n = 1:n1
        coord1[n] = (n-0.5)*Δ1
    end

    coord2 = Array{Float64,1}(undef,n2)
    for n = 1:n2
        coord2[n] = (n-0.5)*Δ2
    end

    coord3 = Array{Float64,1}(undef,n3)
    for n = 1:n3
        coord3[n] = (n-0.5)*Δ3
    end

    if useCUDA
        # Convert to CUDA Arrays
        println(typeof(coord1))
        coord1 = adapt(CuArray,coord1)
        coord2 = adapt(CuArray,coord2)
        coord3 = adapt(CuArray,coord3)
        println(typeof(coord1))
    end
    return coord1, coord2, coord3
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
