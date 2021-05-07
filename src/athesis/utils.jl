using TimerOutputs
using CUDA

macro synctimeit(timer, label, body)
    if has_cuda()
        return :( 
            if ($(esc(timer)).enabled) 
                @timeit $(esc(timer)) $(esc(label)) CUDA.@sync $(esc(body)) 
            else
                $(esc(body))
            end
            )
    else
        return :( 
            if ($(esc(timer)).enabled) 
                @timeit $(esc(timer)) $(esc(label)) $(esc(body)) 
            else
                $(esc(body))
            end
            )
    end

end

macro withCUDA(expr)
    return has_cuda() ? :($(esc(expr))) : :(nothing)
end
