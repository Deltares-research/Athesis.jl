using CUDA

function runSimulation!(simulation)
    runSimulation!(simulation, TimerOutput())
end

function runSimulation!(simulation, to::TimerOutput)

    @synctimeit to "run simulation" begin

        # Unpack time data
        timeData = simulation.timeData
        maxsteps = timeData.maxsteps

        # Now start the time loop
        println("Starting time loop ...")

        kiter = 0
        for n = 1:maxsteps

            simulationConverged, Δh_max, max_index = doExplicitIter!(n, simulation, to::TimerOutput)
            kiter += 1

            if simulationConverged
                println()
                println("-- simulation converged --")
                println("Δh_max = ", Δh_max, " at ", max_index)
                println("nr. of iters = ", kiter)
                println()

                return true, kiter
            end
        end

        return false, kiter
    end # end timer
end

function doTimeStep!(simulation)
    doTimeStep!(simulation, TimerOutput())
end

function doTimeStep!(simulation, to)    
    # Unpack data   
    state       = simulation.state
    timeData    = simulation.timeData

    Δt          = timeData.Δt
    time        = timeData.time

    prepareAndSolve!(simulation, to)

    # Update the old to the new solution
    @synctimeit to "update state" begin
        updateState!(state)
        simulation.timeData.time += Δt
    end

    nothing
end

function doExplicitIter!(n, simulation)
    doExplicitIter!(n, simulation, TimerOutput())
end

function doExplicitIter!(n, simulation, to)

    # Unpack data   
    state       = simulation.state
    timeData    = simulation.timeData
    solverData  = simulation.solverData

    Δt          = timeData.Δt
    time        = timeData.time

    prepareAndSolve!(simulation, to)

    # check convergence
    @synctimeit to "convergence check" begin
        simulationConverged, Δh_max, max_index = checkConvergence!(state, solverData)
    end

    @synctimeit to "print iterations" begin
        if mod(n, 100) == 0
            println("iter = ", n, ", Δh_max = ", Δh_max, " at ", max_index)
        end
    end

    # Update the old to the new solution
    @synctimeit to "update state" begin
        updateState!(state)
        time += Δt
    end

    return simulationConverged, Δh_max, max_index
end

function prepareAndSolve!(simulation, to)
        # Unpack data
        grid        = simulation.grid
        model       = simulation.model
        state       = simulation.state
        parameters  = simulation.parameters
        timeData    = simulation.timeData
        time        = timeData.time
    
        # Add the sources
        @synctimeit to "set rhs" begin
            setSources!(time, model.source)
            setRecharge!(time, model.recharge)
        end
    
        @synctimeit to "set_boundaries" begin
        bc = model.boundaryConditions
            setBoundaries!(grid, state, bc)
        end
    
            @synctimeit to "solve pressure" begin
            pressureEquation!(grid, model, state, parameters, timeData)
        end
    
        @synctimeit to "solve darcy" begin
            darcyEquation!(grid, model, state, parameters, timeData)
        end    
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
    m = reduce(max, h)
    return m, -1
end
