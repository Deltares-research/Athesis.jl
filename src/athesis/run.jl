using CUDA

function runSimulation!(simulation, to::TimerOutput)

    # Unpack time data
    timeData = simulation.timeData
    maxsteps = timeData.maxsteps

    # Now start the time loop
    println("Starting time loop ...")

    for n = 1:maxsteps

        simulationConverged, Δh_max, max_index = doTimestep!(n, simulation, to::TimerOutput)

        if simulationConverged
            println()
            println("-- simulation converged --")
            println("Δh_max = ", Δh_max, " at ", max_index)
            println("nr. of iters = ", n)
            println()
            break
        end
    end
end

function doTimestep!(n, simulation, to::TimerOutput)

    # Unpack data
    grid        = simulation.grid
    model       = simulation.model
    state       = simulation.state
    parameters  = simulation.parameters
    timeData    = simulation.timeData
    solverData  = simulation.solverData

    Δt          = timeData.Δt
    tend        = timeData.tend
    maxsteps    = timeData.maxsteps
    time        = timeData.time

    # Add the sources
    #@timeit to "set rhs" CUDA.@sync begin
    @timeit to "set rhs" begin
        setSources!(time, model.source)
        setRecharge!(time, model.recharge)
    end

    #@timeit to "set_boundaries" CUDA.@sync begin
        @timeit to "set_boundaries" begin
        bc = model.boundaryConditions
        setBoundaries!(grid, state, bc)
    end

    #@timeit to "solve pressure" CUDA.@sync begin
    @timeit to "solve pressure" begin
        pressureEquation!(grid, model, state, parameters, timeData)
    end

    #@timeit to "solve darcy" CUDA.@sync begin
    @timeit to "solve darcy" begin
        darcyEquation!(grid, model, state, parameters, timeData)
    end

    # check convergence
    #@timeit to "convergence check" CUDA.@sync begin
    @timeit to "convergence check" begin
        simulationConverged, Δh_max, max_index = checkConvergence!(state, solverData)
    end

    #@timeit to "print iterations" CUDA.@sync begin
    @timeit to "print iterations" begin
        if mod(n,100) == 0
            println("iter = ", n, ", Δh_max = ", Δh_max, " at ", max_index)
        end
    end

    # Update the old to the new solution
    #@timeit to "update state" CUDA.@sync begin
    @timeit to "update state" begin
        updateState!(state)
        time += Δt
    end

    return simulationConverged, Δh_max, max_index
end

function updateState!(state)
    state.h.parent .= state.hⁿ⁺¹.parent
    state.u.parent .= state.uⁿ⁺¹.parent
    state.v.parent .= state.vⁿ⁺¹.parent
    state.w.parent .= state.wⁿ⁺¹.parent
end

function checkConvergence!(state, solverData)
    solverData.Δh.parent .= abs.(state.hⁿ⁺¹.parent - state.h.parent)
    Δh_max, max_index = findMaximum(solverData.Δh.parent)
    if Δh_max > solverData.ΔhConv
        return false, Δh_max, max_index
    end

    # convergence
    return true, Δh_max, max_index
end

function findMaximum(h)
    return findmax(h)
end

function findMaximum(h::CuArray)
    m = reduce(max,h)
    return m,-1
end
