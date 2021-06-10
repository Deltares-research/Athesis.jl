# initialize.jl

using CUDA
using OffsetArrays
using TimerOutputs

""" Initialize the state vector of the model."""
function initModelState(grid, h0, u0, v0, w0, useCUDA, myFloat)
    # Allocate model state
    # with 4 parameters
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    # Generate initial arrays with start index 0
    useOffset = true

    h    = initField(h0, nx, ny, nz, useCUDA, useOffset, myFloat)
    u    = initField(u0, nx, ny, nz, useCUDA, useOffset, myFloat)
    v    = initField(v0, nx, ny, nz, useCUDA, useOffset, myFloat)
    w    = initField(w0, nx, ny, nz, useCUDA, useOffset, myFloat)

    hⁿ⁺¹ = initField(h0, nx, ny, nz, useCUDA, useOffset, myFloat)
    uⁿ⁺¹ = initField(u0, nx, ny, nz, useCUDA, useOffset, myFloat)
    vⁿ⁺¹ = initField(v0, nx, ny, nz, useCUDA, useOffset, myFloat)
    wⁿ⁺¹ = initField(w0, nx, ny, nz, useCUDA, useOffset, myFloat)

    state = State(h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
    return state
end

""" Instantiate a simulation object (without timer output)."""
function initSimulation(modelInput, useCUDA, myFloat)
    initSimulation(modelInput, useCUDA, myFloat, TimerOutput())
end

""" Instantiate a simulation object (with timer output)."""
function initSimulation(modelInput, useCUDA, myFloat, to)

    @synctimeit to "initialization" begin
        # Initialize the correct sizes/dimensions of the model
        nx        = modelInput.nx
        ny        = modelInput.ny
        nz        = modelInput.nz
        Δx        = modelInput.Δx
        Δy        = modelInput.Δy
        Δz        = modelInput.Δz
        Δt        = modelInput.Δt
        tend      = modelInput.tend
        K0        = modelInput.K0
        S0        = modelInput.S0
        h0        = modelInput.h0
        u0        = modelInput.u0
        v0        = modelInput.v0
        w0        = modelInput.w0
        source    = modelInput.source
        i_src     = modelInput.i_src
        j_src     = modelInput.j_src
        k_src     = modelInput.k_src
        duration  = modelInput.duration
        ΔhConv    = modelInput.ΔhConv
        constRecharge = modelInput.constRecharge
        rechargeFactor = modelInput.rechargeFactor
        boundaryPressure = modelInput.boundaryPressure

        useOffset = true
        noOffset  = false

        # The grid
        x, y, z    = gridCoords(nx, ny, nz, Δx, Δy, Δz, useCUDA, useOffset, myFloat)
        AT         = typeof(x)
        grid       = Grid{AT,myFloat}(nx, ny, nz, Δx, Δy, Δz, x, y, z, useCUDA ? GPUArch() : CPUArch())

        # State vector
        state      = initModelState(grid, h0, u0, v0, w0, useCUDA, myFloat)

        # External forcing
        externals  = initField(0.0, nx, ny, nz, useCUDA, noOffset, myFloat)

        # Model parameters
        K          = initField(K0, nx, ny, nz, useCUDA, useOffset, myFloat)
        specificStorage = S0

        # Sources/sinks
        AT = typeof(externals)
        source     = Source{AT,myFloat}(i_src, j_src, k_src, duration, source, externals)
        setSources!(0.0, source)

        # Recharge
        recharge = Recharge{myFloat}(constRecharge, 0.0, rechargeFactor)

        # Store the boundary conditions
        if (useCUDA)
            boundaryPressure = CuArray(boundaryPressure)
        end
        boundaryConditions = BoundaryConditions(boundaryPressure)


        # Group some parameters in the model.
        # For now sources and boundary conditions
        model      = Model(source, recharge, boundaryConditions)

        # Initialize the set of parameters (for now only K)
        AT = typeof(K)
        parameters = Parameters{AT,myFloat}(K, specificStorage)

        # Input object is no longer needed
        modelInput = nothing

        # Time related data
        maxsteps  = round(Int64, tend / Δt)
        time      = myFloat(0.0)
        timeData  = TimeData(Δt, tend, time, maxsteps)

        # Solver data
        Δh = initField(0.0, nx, ny, nz, useCUDA, useOffset, myFloat)
        AT = typeof(Δh)
        solverData = SolverData{AT,myFloat}(ΔhConv, Δh)

        return Simulation(grid, model, state, parameters, timeData, solverData)
    end # end timer
end
