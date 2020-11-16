# grids.jl

using OffsetArrays

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

function gridCoords(nx, ny, nz, Δx, Δy, Δz, useCUDA, useOffset)
    # Create 3D structured grid
    # For now cell centered

    x = initField1D(0.0, nx, useCUDA, useOffset)
    y = initField1D(0.0, ny, useCUDA, useOffset)
    z = initField1D(0.0, nz, useCUDA, useOffset)

    # Now fill the grid values
    for n = 0:nx+1
        x[n] = (n-0.5)*Δx
    end

    for n = 0:ny+1
        y[n] = (n-0.5)*Δy
    end

    for n = 0:nz+1
        z[n] = (n-0.5)*Δz
    end

    return x, y, z
end
