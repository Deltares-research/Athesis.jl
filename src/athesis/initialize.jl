# initialize.jl

using CuArrays
using Adapt

include("model_input.jl")
include("grids.jl")
include("external_forcing.jl")
include("model.jl")


function init_model_state(grid::Grid2dh, h0, u0, v0, useCUDA)
    # Allocate model state for a 2DH model
    # with 3 parameters
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
    nx = grid.nx
    ny = grid.ny
    h = fill(h0, (nx, ny))
    u = fill(u0, (nx, ny))
    v = fill(v0, (nx, ny))

    if useCUDA
        # Convert to CUDA Arrays
        h = adapt(CuArray,h)
        u = adapt(CuArray,u)
        v = adapt(CuArray,v)
    end

    # Copy the state to the updated state
    hⁿ⁺¹ = copy(h)
    uⁿ⁺¹ = copy(u)
    vⁿ⁺¹ = copy(v)

    state = State2dh(h, u, v, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹)
    return state
end

function init_model_state(grid::Grid2dv, h0, u0, w0, useCUDA)
    # Allocate model state for a 2DV model
    # with 3 parameters
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
    nx = grid.nx
    nz = grid.nz
    h = fill(h0, (nx, nz))
    u = fill(u0, (nx, nz))
    w = fill(w0, (nx, nz))

    if useCUDA
        # Convert to CUDA Arrays
        h = adapt(CuArray,h)
        u = adapt(CuArray,u)
        w = adapt(CuArray,w)
    end

    # Copy the state to the updated state
    hⁿ⁺¹ = copy(h)
    uⁿ⁺¹ = copy(u)
    wⁿ⁺¹ = copy(w)

    state = State2dv(h, u, w, hⁿ⁺¹, uⁿ⁺¹, wⁿ⁺¹)
    return state
end

function init_model_state(grid::Grid3d, h0, u0, v0, w0, useCUDA)
    # Allocate model state for a 3D model
    # with 4 parameters
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
     nx = grid.nx
    ny = grid.ny
    nz = grid.nz
    h = fill(h0, (nx, ny, nz))
    u = fill(u0, (nx, ny, nz))
    v = fill(v0, (nx, ny, nz))
    w = fill(w0, (nx, ny, nz))

    if useCUDA
        # Convert to CUDA Arrays
        h = adapt(CuArray,h)
        u = adapt(CuArray,u)
        v = adapt(CuArray,v)
        w = adapt(CuArray,w)
    end

    # Copy the state to the updated state
    hⁿ⁺¹ = copy(h)
    uⁿ⁺¹ = copy(u)
    vⁿ⁺¹ = copy(v)
    wⁿ⁺¹ = copy(w)

    state = State3d(h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
    return state
end

function init_parameters(n1, n2, p0, useCUDA)
    # Allocate parameters for a 2D model (2DH or 2DV)
    # with 1 parameter
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
    p1 = fill(p0, (n1, n2))

    if useCUDA
        # Convert to CUDA Array
        p1 = adapt(CuArray,p1)
    end
    return p1
end

function init_parameters(n1, n2, n3, p0, useCUDA)
    # Allocate parameters for a 3D model
    # with 1 parameter
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
    p1 = fill(p0, (n1, n2, n3))
    if useCUDA
        # Convert to CUDA Array
        p1 = adapt(CuArray,p1)
    end

    return p1
end



function model_initialize()

    println("Initialize the model ...")

    # Get/read model input
    input = model_input()

    # Set model type/dimension
    model_type = input.model_type

    # Initialize the correct sizes/dimensions of the model
    if model_type == "2DH"

        #Data_type = input.data_type
        useCUDA   = input.useCUDA
        ndim      = input.ndim
        nx        = input.nx
        ny        = input.ny
        Δx        = input.Δx
        Δy        = input.Δy
        Δt        = input.Δt
        tend      = input.tend
        K0        = input.K0
        h0        = input.h0
        u0        = input.u0
        v0        = input.v0
        source    = input.source
        i_src     = input.i_src
        j_src     = input.j_src
        duration  = input.duration

        # The grid
        x, y       = grid_coords(nx,ny, Δx, Δy, useCUDA)
        grid       = Grid2dh(ndim, nx, ny, Δx, Δy, x, y)

        # Model state vector
        state      = init_model_state(grid, h0, u0, v0, useCUDA)

        # External forcing
        externals  = init_externals(nx, ny, useCUDA)

        # Model parameters
        K          = init_parameters(nx, ny, K0, useCUDA)

        # Sources/sinks
        source     = Source2dh(i_src, j_src, duration, source, externals)
        set_sources!(0.0, source)

    elseif model_type == "2DV"

        #Data_type = input.Data_type
        useCUDA   = input.useCUDA
        ndim      = input.ndim
        nx        = input.nx
        nz        = input.nz
        Δx        = input.Δx
        Δz        = input.Δz
        Δt        = input.Δt
        tend      = input.tend
        K0        = input.K0
        h0        = input.h0
        u0        = input.u0
        w0        = input.w0
        source    = input.source
        i_src     = input.i_src
        k_src     = input.k_src
        duration  = input.duration

        # The grid
        x, z       = grid_coords(nx, nz, Δx, Δz, useCUDA)
        grid       = Grid2dv(ndim, nx, nz, Δx, Δz, x, z)

        # Model state vector
        state      = init_model_state(grid, h0, u0, w0, useCUDA)

        # External forcing
        externals  = init_externals(nx, nz, useCUDA)

        # Model parameters
        K          = init_parameters(nx, nz, K0, useCUDA)

        # Sources/sinks
        source     = Source2dv(i_src, k_src, duration, source, externals)
        set_sources!(0.0, source)

    elseif model_type == "3D"

        #Data_type = input.Data_type
        useCUDA   = input.useCUDA
        ndim      = input.ndim
        nx        = input.nx
        ny        = input.ny
        nz        = input.nz
        Δx        = input.Δx
        Δy        = input.Δy
        Δz        = input.Δz
        Δt        = input.Δt
        tend      = input.tend
        K0        = input.K0
        h0        = input.h0
        u0        = input.u0
        v0        = input.v0
        w0        = input.w0
        source    = input.source
        i_src     = input.i_src
        j_src     = input.j_src
        k_src     = input.k_src
        duration  = input.duration

        # The grid
        x, y, z    = grid_coords(nx, ny, nz, Δx, Δy, Δz, useCUDA)
        grid       = Grid3d(ndim, nx, ny, nz, Δx, Δy, Δz, x, y, z)

        # State vector
        state      = init_model_state(grid, h0, u0, v0, w0, useCUDA)

        # External forcing
        externals  = init_externals(nx, ny, nz, useCUDA)

        # Model parameters
        K          = init_parameters(nx, ny, nz, K0, useCUDA)

        # Sources/sinks
        source     = Source3d(i_src, j_src, k_src, duration, source, externals)
        set_sources!(0.0, source)

    else
        throw("ERROR: erroneous modeltype: ", model_type, ". Allowed model_types are:\n\t2DH\n\t2DV\n\t3D")
    end

    # Input object is no longer needed
    input = nothing

    # Group some parameters in the model.
    # For now only sources
    model      = Model(model_type, source)

    # Initialize the set of parameters (for now only K)
    parameters = Parameters(K)

    # Time related data
    maxsteps   = round(Int64, tend/Δt)
    time       = 0.0
    time_data  = Time_data(Δt, tend, time, maxsteps)

    return grid, model, state, parameters, time_data
end
