# fields.jl

using CUDA
using OffsetArrays
using Adapt

Adapt.adapt_structure(to, x::OffsetArray) = OffsetArray(adapt(to, parent(x)), x.offsets)


function initField(inival, nx, ny, nz, useCUDA, useOffset, myFloat)
    # Initialize a field/array

    # Convert inival to correct Float type
    inival = myFloat(inival)

    # If it is an offset array, increase size by 2
    if useOffset
        var = fill(inival, (nx + 2, ny + 2, nz + 2))
    else
        var = fill(inival, (nx, ny, nz))
    end

    # CUDA array y/n
    if useCUDA
        # Convert to CUDA Arrays
        var = CuArray(var)
    end

    # When offset array, shift start index to 0
    if useOffset
        var = OffsetArray(var, (0:nx + 1, 0:ny + 1, 0:nz + 1))
    end

    return var
end

# fields.jl

using CUDA
using OffsetArrays

function initField1D(inival, nx, useCUDA, useOffset, myFloat)
    # Initialize a field/array

    # Convert inival to correct Float type
    inival = myFloat(inival)

    # If it is an offset array, increase size by 2
    if useOffset
        var = fill(inival, nx + 2)
    else
        var = fill(inival, nx)
    end

    # CUDA array y/n
    if useCUDA
        # Convert to CUDA Arrays
        var = CuArray(var)
    end

    # When offset array, shift start index to 0
    if useOffset
        var = OffsetArray(var, 0:nx + 1)
    end

    return var
end
