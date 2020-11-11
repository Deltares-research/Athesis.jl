using CUDA

function runSimulation(simulation, to::TimerOutput)

    # Unpack data
    grid = simulation.grid
    model = simulation.model
    state = simulation.state
    parameters = simulation.parameters
    time_data = simulation.timeData
    solver_data = simulation.solverData

    Δt       = time_data.Δt
    tend     = time_data.tend
    maxsteps = time_data.maxsteps
    time     = time_data.time

    # Now start the time loop
    println("Starting time loop ...")

    # Initialize the present time
    time = 0.0

    for n = 1:maxsteps

        # Add the sources
        @timeit to "set rhs" CUDA.@sync begin
            set_sources!(time, model.source)
            set_recharge!(time,model.recharge)
        end

        @timeit to "set_boundaries" CUDA.@sync begin
            bc = model.boundary_conditions
            set_boundaries!(grid, state, bc)
        end

        @timeit to "solve pressure" CUDA.@sync begin
            pressure_equation!(grid, model, state, parameters, time_data)
        end
        @timeit to "solve darcy" CUDA.@sync begin
            darcy_equation!(grid, model, state, parameters, time_data)
        end

        # check convergence
        @timeit to "convergence check" CUDA.@sync begin
            has_converged, Δh_max, max_index = check_convergence!(state, solver_data)
            if (has_converged)
                println()
                println("-- solver converged --")
                println("Δh_max = ", Δh_max, " at ", max_index)
                println("nr. of iters = ", n)
                println()

                # Update the old to the new solution
                updateState!(state)
                time += Δt

                break
            end
        end

        @timeit to "print iterations" CUDA.@sync begin
            if mod(n,100) == 0
                println("iter = ", n, ", Δh_max = ", Δh_max, " at ", max_index)
            end
        end

        # Update the old to the new solution
        @timeit to "update state" CUDA.@sync begin
            updateState!(state)
            time += Δt
        end

    end
end

function updateState!(state)
    state.h.parent .= state.hⁿ⁺¹.parent
    state.u.parent .= state.uⁿ⁺¹.parent
    state.v.parent .= state.vⁿ⁺¹.parent
    state.w.parent .= state.wⁿ⁺¹.parent
end

function check_convergence!(state, solver_data)
    solver_data.Δh.parent .= abs.(state.hⁿ⁺¹.parent - state.h.parent)
    Δh_max, max_index = find_maximum(solver_data.Δh.parent)
    if Δh_max > solver_data.hclose
        return false, Δh_max, max_index
    end

    # convergence
    return true, Δh_max, max_index
end

function find_maximum(h)
    return findmax(h)
end

function find_maximum(h::CuArray)
    m = reduce(max,h)
    return m,-1
end
