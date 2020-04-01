using CUDAdrv, CUDAnative
using CuArrays

function cuda_wrap_kernel!(kernel, values::AbstractArray{Float64}, result::AbstractArray{Float64})
    ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
    if (ix <= length(values))
        @inbounds result[ix] = kernel(values, ix)
    end

    return nothing
end

function gridloop(kernel, field, result)
    if field.values isa CuArray
        ths = 256
        bls = Int(ceil(field.grid.N / ths))
        @cuda threads=ths blocks=bls cuda_wrap_kernel!(kernel, field.values, result)
    else
        for i = 1 : field.grid.N
            @inbounds result[i] = kernel(field.values, i)
        end
    end
end

function run()

    println("Hit c for cuda...")
    c = readline()
    useCUDA = c == "c"

    grid = CartesianGrid(10, 100.0, 0.0)

    fieldValues = collect(1.0:grid.N)
    if useCUDA
        fieldValues = CuArray(fieldValues)
    end
    field = Field(fieldValues, grid)

    result = similar(field.values)

    println("before: ", field.values)
    gridloop(kernel, field, result)
    println("after:  ", result)
end

mutable struct CartesianGrid
    N :: Int64   # number of cells
    L :: Float64 # domain size
    Δ :: Float64 # grid spacing
end

function CartesianGrid(nCells, length)
    return CartesianGrid(nCells, length, length/nCells)
end

mutable struct Field
    values :: AbstractArray{Float64}
    grid :: CartesianGrid
end


@inline function kernel(ξ, i)
    return 2 * ξ[i]
end
