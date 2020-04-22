# initialize.jl

using CuArrays
using Adapt

include("model_input.jl")
include("grids.jl")
include("external_forcing.jl")
include("model.jl")


function init_model_state(grid::Grid2dh, h0, u0, v0, Data_type)
    # Allocate model state for a 2DH model
    # with 3 parameters
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
    nx = grid.nx
    ny = grid.ny
    h = fill(h0, (nx, ny))
    u = fill(u0, (nx, ny))
    v = fill(v0, (nx, ny))

    # Convert to requested data_type (possibly CuArray)
    h = adapt(Data_type,h)
    u = adapt(Data_type,u)
    v = adapt(Data_type,v)

    # Copy the state to the updated state
    hⁿ⁺¹ = copy(h)
    uⁿ⁺¹ = copy(u)
    vⁿ⁺¹ = copy(v)

    #state = State2dh{data_type}(h, u, v, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹)
    state = State2dh(h, u, v, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹)
    return state
end

function init_model_state(grid::Grid2dv, h0, u0, w0, Data_type)
    # Allocate model state for a 2DV model
    # with 3 parameters
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
    nx = grid.nx
    nz = grid.nz
    h = fill(h0, (nx, nz))
    u = fill(u0, (nx, nz))
    w = fill(w0, (nx, nz))

    # Convert to requested data_type (possibly CuArray)
    h = adapt(Data_type,h)
    u = adapt(Data_type,u)
    w = adapt(Data_type,w)

    # Copy the state to the updated state
    hⁿ⁺¹ = copy(h)
    uⁿ⁺¹ = copy(u)
    wⁿ⁺¹ = copy(w)

    #state = State2dv{data_type}(h, u, w, hⁿ⁺¹, uⁿ⁺¹, wⁿ⁺¹)
    state = State2dv(h, u, w, hⁿ⁺¹, uⁿ⁺¹, wⁿ⁺¹)
    return state
end

function init_model_state(grid::Grid3d, h0, u0, v0, w0, Data_type)
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

    # Convert to requested data_type (possibly CuArray)
    h = adapt(Data_type,h)
    u = adapt(Data_type,u)
    v = adapt(Data_type,v)
    w = adapt(Data_type,w)

    # Copy the state to the updated state
    hⁿ⁺¹ = copy(h)
    uⁿ⁺¹ = copy(u)
    vⁿ⁺¹ = copy(v)
    wⁿ⁺¹ = copy(w)

    #state = State3d{data_type}(h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
    state = State3d(h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
    return state
end

function init_parameters(n1, n2, p0, Data_type)
    # Allocate parameters for a 2D model (2DH or 2DV)
    # with 1 parameter
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
    p1 = fill(p0, (n1, n2))
    p1 = adapt(Data_type,p1)
    return p1
end

function init_parameters(n1, n2, n3, p0, Data_type)
    # Allocate parameters for a 3D model
    # with 1 parameter
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
    p1 = fill(p0, (n1, n2, n3))
    p1 = adapt(Data_type,p1)
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

        Data_type = input.data_type
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
        x, y       = grid_coords(nx,ny, Δx, Δy, Data_type)
        #grid       = Grid2dh{data_type}(ndim, nx, ny, Δx, Δy, x, y)
        grid       = Grid2dh(ndim, nx, ny, Δx, Δy, x, y)

        # Model state vector
        state      = init_model_state(grid, h0, u0, v0, Data_type)

        # External forcing
        externals  = init_externals(nx, ny, Data_type)

        # Model parameters
        K          = init_parameters(nx, ny, K0, Data_type)

        # Sources/sinks
        #source     = Source2dh{data_type}(i_src, j_src, n_src, source, externals)
        source     = Source2dh(i_src, j_src, duration, source, externals)
        set_sources!(0.0, source)

    elseif model_type == "2DV"

        Data_type = input.Data_type
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
        x, z       = grid_coords(nx, nz, Δx, Δz, Data_type)
        #grid       = Grid2dv{data_type}(ndim, nx, nz, Δx, Δz, x, z)
        grid       = Grid2dv(ndim, nx, nz, Δx, Δz, x, z)

        # Model state vector
        state      = init_model_state(grid, h0, u0, w0, Data_type)

        # External forcing
        externals  = init_externals(nx, nz, Data_type)

        # Model parameters
        K          = init_parameters(nx, nz, K0, Data_type)

        # Sources/sinks
        #source     = Source2dv{3}(i_src, k_src, n_src, source, externals)
        source     = Source2dv(i_src, k_src, duration, source, externals)
        set_sources!(0.0, source)

    elseif model_type == "3D"

        Data_type = input.Data_type
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
        x, y, z    = grid_coords(nx, ny, nz, Δx, Δy, Δz, Data_type)
        #grid       = Grid3d{data_type}(ndim, nx, ny, nz, Δx, Δy, Δz, x, y, z)
        grid       = Grid3d(ndim, nx, ny, nz, Δx, Δy, Δz, x, y, z)

        # State vector
        state      = init_model_state(grid, h0, u0, v0, w0, Data_type)

        # External forcing
        externals  = init_externals(nx, ny, nz, Data_type)

        # Model parameters
        K          = init_parameters(nx, ny, nz, K0, Data_type)

        # Sources/sinks
        #source     = Source3d{data_type}(i_src, j_src, k_src, n_src, source, externals)
        source     = Source3d(i_src, j_src, k_src, duration, source, externals)
        set_sources!(0.0, source)

    else
        throw("ERROR: erroneous modeltype: ", model_type, ". Allowed model_types are:\n\t2DH\n\t2DV\n\t3D")
    end

    #
    # Group some parameters in the model.
    # For now only sources
    #model      = Model{data_type}(model_type, source)
    model      = Model(model_type, source)

    #parameters = Parameters{data_type}(K)
    parameters = Parameters(K)

    # Time related data
    maxsteps   = round(Int64, tend/Δt)
    time       = 0.0
    time_data  = Time_data(Δt, tend, time, maxsteps)

    return grid, model, state, parameters, time_data
end
