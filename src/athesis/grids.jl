# grids.jl

using OffsetArrays
using CUDA
using KernelAbstractions
using CUDAKernels

abstract type AbstractArchitecture end
struct CPUArch <: AbstractArchitecture end
struct GPUArch <: AbstractArchitecture end

device(::CPUArch) = KernelAbstractions.CPU()
device(::GPUArch) = CUDAKernels.CUDADevice()

mutable struct Grid{AT,FT}
    # grid_type::String
    nx::Int64
    ny::Int64
    nz::Int64
    Δx::FT
    Δy::FT
    Δz::FT
    x::AT
    y::AT
    z::AT
    arch::AbstractArchitecture
end

function gridCoords(nx, ny, nz, Δx, Δy, Δz, useCUDA, useOffset, myFloat)
    # Create 3D structured grid
    # For now cell centered

    x = initField1D(0.0, nx, useCUDA, useOffset, myFloat)
    y = initField1D(0.0, ny, useCUDA, useOffset, myFloat)
    z = initField1D(0.0, nz, useCUDA, useOffset, myFloat)

    # Now fill the grid values
    for n = 0:nx + 1
        x[n] = myFloat(n - 0.5) * Δx
    end

    for n = 0:ny + 1
        y[n] = myFloat(n - 0.5) * Δy
    end

    for n = 0:nz + 1
        z[n] = myFloat(n - 0.5) * Δz
    end

    return x, y, z
end
