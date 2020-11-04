include("model.jl")

function init_externals(nx, ny, nz, useCUDA)
    externals = zeros(nx, ny, nz)
    if useCUDA
        # Convert to CUDA Array
        externals = CuArray(externals)
    end
    return externals
end


function set_sources!(time, source)
    # 3D implementation
    i_src = source.i_src
    j_src = source.j_src
    k_src = source.k_src
    duration = source.duration
    source.external_source .= 0.0
    if time < duration
        # Set a single source of at the prescribed point
        source.external_source[i_src, j_src, k_src] = source.discharge
    else
        source.external_source[i_src, j_src, k_src] = 0.0
    end
end


function set_recharge!(time, recharge)
    #duration = recharge.duration
    #if time < duration
        # Set a single source of at the prescribed point
        recharge.recharge_flux = recharge.recharge_factor*recharge.const_recharge
    #else
    #    recharge.recharge_flux = 0.0
    #end
end
