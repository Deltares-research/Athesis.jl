using CUDAdrv, CUDAnative
using CuArrays
using TimerOutputs

function backendname(useCUDA)
    return useCUDA ? "GPU" : "CPU"
end

function run()

    to = TimerOutput()

    N = 100000000

    println("Hit c for cuda...")
    c = readline()
    useCUDA = (c == "c")

    fieldA = fill(1.0, N)
    fieldB = fill(2.0, N)
    if useCUDA
        fieldA = CuArray(fieldA)
        fieldB = CuArray(fieldB)
    end
    result = similar(fieldA)

    # first separate kernels:
    fill!(result, 0.0)
    @timeit to "separate "*backendname(useCUDA) begin
    gridloop(result, kernel_multby2, fieldA)
    gridloop(result, kernel_mean, (fieldA, fieldB))
    end
    println("result separate:  ", result[1:10])

    # then merged into one:
    fill!(result, 0.0)
    @timeit to "fused "*backendname(useCUDA) begin
    gridloop(result, kernel_fused, (fieldA, fieldB))
    end
    println("result fused:  ", result[1:10])

    # then composed:
    fill!(result, 0.0)
    @timeit to "composed "*backendname(useCUDA) begin
    gridloop(result, kernel_composed, (fieldA, fieldB))
    end
    println("result composed:  ", result[1:10])

    print_timer(to)

end

@inline function kernel_multby2(i, result, ξ)
    @inbounds result[i] += 2 * ξ[i]
end

@inline function kernel_mean(i, result, (ξ, ψ))
    @inbounds result[i] += 0.5 * (ξ[i] + ψ[i])
end

@inline function kernel_fused(i, result, (ξ, ψ))
    @inbounds result[i] += 2*ξ[i] + 0.5*(ξ[i] + ψ[i])
end

@inline function kernel_composed(i, result, (ξ, ψ))
    kernel_multby2(i, result, ξ)
    kernel_mean(i, result, (ξ, ψ))
end

function cuda_wrap_kernel(result, kernel, input)
    ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
    if (ix <= length(result))
        kernel(ix, result, input)
    end

    return nothing
end

function gridloop(result, kernel, input)
    if result isa CuArray
        ths = 256
        bls = Int(ceil(length(result) / ths))
        @cuda threads=ths blocks=bls cuda_wrap_kernel(result, kernel, input)
    else
        for i = 1 : length(result)
            kernel(i, result, input)
        end
    end
end
