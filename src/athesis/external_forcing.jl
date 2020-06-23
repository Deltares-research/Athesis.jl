# external_forcing.jl

using Adapt

include("model.jl")

function init_externals(n1, n2, useCUDA)
    externals = zeros(n1, n2)
    if useCUDA
        # Convert to CUDA Array
        externals = adapt(CuArray,externals)
    end

    return externals
end

function init_externals(n1, n2, n3, useCUDA)
    externals = zeros(n1, n2, n3)
    if useCUDA
        # Convert to CUDA Array
        externals = adapt(CuArray,externals)
    end
    return externals
end

function set_sources!(time::Float64, source::Source2dh)
    # 2DH implementation
    i_src = source.i_src
    j_src = source.j_src
    duration = source.duration
    if time < duration
        # Set a single source at the prescribed point
        source.external_source[i_src,j_src] = source.discharge
    else
        source.external_source[i_src,j_src] = 0.0
    end
end

function set_sources!(time::Float64, source::Source2dv)
    # 2DV implementation
    i_src = source.i_src
    k_src = source.k_src
    duration = source.duration
    if time < duration
        # Set a single source at the prescribed point
        source.external_source[i_src,k_src] = source.discharge
    else
        source.external_source[i_src,k_src] = 0.0
    end
end

function set_sources!(time::Float64, source::Source3d)
    # 3D implementation
    i_src = source.i_src
    j_src = source.j_src
    k_src = source.k_src
    duration = source.duration
    if time < duration
        # Set a single source of at the prescribed point
        source.external_source[i_src, j_src, k_src] = source.discharge
    else
        source.external_source[i_src, j_src, k_src] = 0.0
    end
end
