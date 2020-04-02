using CUDAdrv, CUDAnative
using CuArrays
using TimerOutputs

function backendname(useCUDA)
    return useCUDA ? "GPU" : "CPU"
end

function run()

    to = TimerOutput()

    N = 10000000

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
    println("separate, in: ", fieldA[1:10])
    @timeit to "separate "*backendname(useCUDA) begin
    gridloop(kernel_multby2, result, fieldA)
    gridloop(kernel_mean, result, (fieldA, fieldB))
    end
    println("result:  ", result[1:10])

    # then merged into one:
    fill!(result, 0.0)
    println("fused, in: ", fieldA[1:10])
    @timeit to "fused "*backendname(useCUDA) begin
    gridloop(kernel_fused, result, (fieldA, fieldB))
    end
    println("result:  ", result[1:10])

    # then composed:
    fill!(result, 0.0)
    println("composed, in: ", fieldA[1:10])
    @timeit to "composed "*backendname(useCUDA) begin
    gridloop(kernel_composed, result, (fieldA, fieldB))
    end
    println("result:  ", result[1:10])

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

function cuda_wrap_kernel(kernel::Function, result::AbstractArray{Float64}, args)
    ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
    if (ix <= length(result))
        kernel(ix, result, args)
    end

    return nothing
end

function gridloop(kernel::Function, result::AbstractArray{Float64}, args)
    if result isa CuArray
        ths = 256
        bls = Int(ceil(length(result) / ths))
        @cuda threads=ths blocks=bls cuda_wrap_kernel(kernel, result, args)
    else
        for i = 1 : length(result)
            kernel(i, result, args)
        end
    end
end
