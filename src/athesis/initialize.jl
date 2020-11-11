# initialize.jl

using CUDA
using OffsetArrays

function init_model_state(grid, h0, u0, v0, w0, useCUDA)
    # Allocate model state
    # with 4 parameters
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    # Generate initial arrays with start index 1
    h = fill(h0, (nx+2, ny+2, nz+2))
    u = fill(u0, (nx+2, ny+2, nz+2))
    v = fill(v0, (nx+2, ny+2, nz+2))
    w = fill(w0, (nx+2, ny+2, nz+2))

    hn = fill(h0, (nx+2, ny+2, nz+2))
    un = fill(u0, (nx+2, ny+2, nz+2))
    vn = fill(v0, (nx+2, ny+2, nz+2))
    wn = fill(w0, (nx+2, ny+2, nz+2))

    if useCUDA
        # Convert to CUDA Arrays
        h = CuArray(h)
        u = CuArray(u)
        v = CuArray(v)
        w = CuArray(w)

        hn = CuArray(hn)
        un = CuArray(un)
        vn = CuArray(vn)
        wn = CuArray(wn)
    end

    # Shift indices to let the arrays run from index 0
    h = OffsetArray(h, (0:nx+1, 0:ny+1, 0:nz+1))
    u = OffsetArray(u, (0:nx+1, 0:ny+1, 0:nz+1))
    v = OffsetArray(v, (0:nx+1, 0:ny+1, 0:nz+1))
    w = OffsetArray(w, (0:nx+1, 0:ny+1, 0:nz+1))

    hⁿ⁺¹ = OffsetArray(hn, (0:nx+1, 0:ny+1, 0:nz+1))
    uⁿ⁺¹ = OffsetArray(un, (0:nx+1, 0:ny+1, 0:nz+1))
    vⁿ⁺¹ = OffsetArray(vn, (0:nx+1, 0:ny+1, 0:nz+1))
    wⁿ⁺¹ = OffsetArray(wn, (0:nx+1, 0:ny+1, 0:nz+1))

    # Copy the state to the updated state
    #hⁿ⁺¹ = copy(h)
    #uⁿ⁺¹ = copy(u)
    #vⁿ⁺¹ = copy(v)
    #wⁿ⁺¹ = copy(w)

    state = State(h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
    return state
end


function init_parameters(nx, ny, nz, p0, useCUDA)
    # Allocate parameters for a 3D model
    # with 1 parameter
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom

    # Generate initial array with start index 1
    p1 = fill(p0, (nx+2, ny+2, nz+2))

    if useCUDA
        # Convert to CUDA Array
        p1 = CuArray(p1)
    end

    # Shift index to let array start at 0
    p1 = OffsetArray(p1, (0:nx+1, 0:ny+1, 0:nz+1))



    return p1
end


function model_initialize(useCUDA)

    println("Initialize the model ...")

    # Get/read model input
    input = model_input()

    # Initialize the correct sizes/dimensions of the model
    nx        = input.nx
    ny        = input.ny
    nz        = input.nz
    Δx        = input.Δx
    Δy        = input.Δy
    Δz        = input.Δz
    Δt        = input.Δt
    tend      = input.tend
    K0        = input.K0
    S0        = input.S0
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
    boundary_pressure = input.boundary_pressure

    # The grid
    x, y, z    = grid_coords(nx, ny, nz, Δx, Δy, Δz, useCUDA)
    grid       = Grid(nx, ny, nz, Δx, Δy, Δz, x, y, z)

    # State vector
    state      = init_model_state(grid, h0, u0, v0, w0, useCUDA)

    # External forcing
    externals  = init_externals(nx, ny, nz, useCUDA)

    # Model parameters
    K          = init_parameters(nx, ny, nz, K0, useCUDA)
    specific_storage = S0

    # Sources/sinks
    source     = Source(i_src, j_src, k_src, duration, source, externals)
    set_sources!(0.0, source)

    # Recharge
    recharge = Recharge(const_recharge, 0.0, recharge_factor)

    # Store the boundary conditions
    if (useCUDA)
        boundary_pressure = CuArray(boundary_pressure)
    end
    boundary_conditions = BoundaryConditions(boundary_pressure)


    # Group some parameters in the model.
    # For now sources and boundary conditions
    model      = Model(source, recharge, boundary_conditions)

    # Initialize the set of parameters (for now only K)
    parameters = Parameters(K, specific_storage)

    # Input object is no longer needed
    input = nothing

    # Time related data
    maxsteps   = round(Int64, tend/Δt)
    time       = 0.0
    time_data  = Time_data(Δt, tend, time, maxsteps)

    # Solver data
    hclose = 1e-5
    Δh = fill(0.0, (grid.nx+2, grid.ny+2, grid.nz+2))
    if useCUDA
        # Convert to CUDA Arrays
        Δh = CuArray(Δh)
    end
    Δh = OffsetArray(Δh, (0:grid.nx+1, 0:grid.ny+1, 0:grid.nz+1))
    solver_data = Solver_data(hclose, Δh)

    return grid, model, state, parameters, time_data, solver_data
end
