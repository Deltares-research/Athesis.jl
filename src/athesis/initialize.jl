# initialize.jl

using CUDA
using OffsetArrays
using TimerOutputs

function initModelState(grid, h0, u0, v0, w0, useCUDA)
    # Allocate model state
    # with 4 parameters
    # for now cell-centered,
    # i.e. in all direction the same number of degrees of freedom
    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    # Generate initial arrays with start index 0
    useOffset = true

    h    = initField(h0, nx, ny, nz, useCUDA, useOffset)
    u    = initField(u0, nx, ny, nz, useCUDA, useOffset)
    v    = initField(v0, nx, ny, nz, useCUDA, useOffset)
    w    = initField(w0, nx, ny, nz, useCUDA, useOffset)

    hⁿ⁺¹ = initField(h0, nx, ny, nz, useCUDA, useOffset)
    uⁿ⁺¹ = initField(u0, nx, ny, nz, useCUDA, useOffset)
    vⁿ⁺¹ = initField(v0, nx, ny, nz, useCUDA, useOffset)
    wⁿ⁺¹ = initField(w0, nx, ny, nz, useCUDA, useOffset)

    state = State(h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
    return state
end


function initSimulation(modelInput, useCUDA, to::TimerOutput)

    @synctimeit to "initialization" begin

        println("Initialize the model ...")

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
        x, y, z    = gridCoords(nx, ny, nz, Δx, Δy, Δz, useCUDA, useOffset)
        grid       = Grid(nx, ny, nz, Δx, Δy, Δz, x, y, z)

        # State vector
        state      = initModelState(grid, h0, u0, v0, w0, useCUDA)

        # External forcing
        externals  = initField(0.0, nx, ny, nz, useCUDA, noOffset)

        # Model parameters
        K          = initField(K0, nx, ny, nz, useCUDA, useOffset)
        specificStorage = S0

        # Sources/sinks
        source     = Source(i_src, j_src, k_src, duration, source, externals)
        setSources!(0.0, source)

        # Recharge
        recharge = Recharge(constRecharge, 0.0, rechargeFactor)

        # Store the boundary conditions
        if (useCUDA)
            boundaryPressure = CuArray(boundaryPressure)
        end
        boundaryConditions = BoundaryConditions(boundaryPressure)


        # Group some parameters in the model.
        # For now sources and boundary conditions
        model      = Model(source, recharge, boundaryConditions)

        # Initialize the set of parameters (for now only K)
        parameters = Parameters(K, specificStorage)

        # Input object is no longer needed
        modelInput = nothing

        # Time related data
        maxsteps   = round(Int64, tend/Δt)
        time       = 0.0
        timeData  = TimeData(Δt, tend, time, maxsteps)

        # Solver data
        Δh = initField(0.0, nx, ny, nz, useCUDA, useOffset)
        solverData = SolverData(ΔhConv, Δh)

        return Simulation(grid, model, state, parameters, timeData, solverData)
    end # end timer
end
