using TimerOutputs
using CUDA

macro synctimeit(timer, label, body)
    if has_cuda()
        return :( @timeit $(esc(timer)) $(esc(label)) CUDA.@sync $(esc(body)) )
    else
        return :( @timeit $(esc(timer)) $(esc(label)) $(esc(body)) )
    end
end

macro withCUDA(expr)
    return has_cuda() ? :($(esc(expr))) : :(nothing)
end
