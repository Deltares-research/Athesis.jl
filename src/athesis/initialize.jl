# initialize.jl

using CuArrays
using Adapt

include("model_input.jl")
include("grids.jl")
include("external_forcing.jl")
include("model.jl")

function init_model_state(grid, h0, u0, v0, w0, useCUDA)
    # Allocate model state
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
        # h = adapt(CuArray,h)
        # u = adapt(CuArray,u)
        # v = adapt(CuArray,v)
        # w = adapt(CuArray,w)
        h = CuArray(h)
        u = CuArray(u)
        v = CuArray(v)
        w = CuArray(w)
    end

    # Copy the state to the updated state
    hⁿ⁺¹ = copy(h)
    uⁿ⁺¹ = copy(u)
    vⁿ⁺¹ = copy(v)
    wⁿ⁺¹ = copy(w)

    state = State(h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
    return state
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

    # Initialize the correct sizes/dimensions of the model
    useCUDA   = input.useCUDA
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
    const_recharge = input.const_recharge
    recharge_factor = input.recharge_factor

    # The grid
    x, y, z    = grid_coords(nx, ny, nz, Δx, Δy, Δz, useCUDA)
    grid       = Grid(nx, ny, nz, Δx, Δy, Δz, x, y, z)

    # State vector
    state      = init_model_state(grid, h0, u0, v0, w0, useCUDA)

    # External forcing
    externals  = init_externals(nx, ny, nz, useCUDA)

    # Model parameters
    K          = init_parameters(nx, ny, nz, K0, useCUDA)

    # Sources/sinks
    source     = Source(i_src, j_src, k_src, duration, source, externals)
    set_sources!(0.0, source)

    # Recharge
    recharge = Recharge(const_recharge, 0.0, recharge_factor)
    
    # Input object is no longer needed
    input = nothing

    # Group some parameters in the model.
    # For now only sources
    model      = Model(source, recharge)

    # Initialize the set of parameters (for now only K)
    parameters = Parameters(K)

    # Time related data
    maxsteps   = round(Int64, tend/Δt)
    time       = 0.0
    time_data  = Time_data(Δt, tend, time, maxsteps)

    return grid, model, state, parameters, time_data
end
