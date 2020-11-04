# 3D groundwater
using Plots
using CUDA
using TimerOutputs

using Athesis

# This is the present data storage:
# grid       = (nx, ny, nz, Δx, Δy, Δz, x, y, z)
# state      = (h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
# model      = externals, recharge
# parameters = (K, Δt, tend)
# externals  = source
# source     = (i_src, j_src, k_src, n_src, externals)


function groundwater3d(isBenchmark = false, useGPU = false)

    # Backend selection
    useCUDA = false
    if isBenchmark
        useCUDA = useGPU
    else
        println("Hit c for cuda...")
        c = readline()
        useCUDA = (c == "c")
    end

    to = TimerOutput()

    @timeit to "run time" begin

    # Initialize the model
    println("Running 3D groundwater model:")

    @timeit to "initialization" begin
        grid, model, state, parameters, time_data, solver_data = model_initialize(useCUDA)

        # Unpack time parameters required for the time loop
        Δt       = time_data.Δt
        tend     = time_data.tend
        maxsteps = time_data.maxsteps
        time     = time_data.time

        # Now start the time loop
        println("Starting time loop ...")

        # Initialize the present time
        time = 0.0
    end

    @timeit to "solve" begin

    bulge = []

    for n = 1:maxsteps

        # Add the sources
        @timeit to "set rhs" begin
            set_sources!(time, model.source)
            set_recharge!(time,model.recharge)
        end

        @timeit to "set_boundaries" begin
            bc = model.boundary_conditions
            set_boundaries!(grid, state, bc)
        end

        @timeit to "solve pressure" begin
            pressure_equation!(grid, model, state, parameters, time_data)
        end
        @timeit to "solve darcy" begin
            #darcy_equation!(grid, model, state, parameters, time_data)
        end

        # check convergence
        @timeit to "convergence check" begin
            has_converged, Δh_max, max_index = check_convergence!(state, solver_data)
            if (has_converged)
                println()
                println("-- solver converged --")
                println("Δh_max = ", Δh_max, " at ", max_index)
                println("nr. of iters = ", n)
                println()

                # todo: update h_n+1 here
                # ...
                break
            end
        end

        @timeit to "append time history" begin
            append!(bulge, state.hⁿ⁺¹[Int(ceil(grid.nx/2)),Int(ceil(grid.ny/2)),grid.nz])
        end

        #if mod(n,100) == 0
        #    plot_model(grid, state)
        #end
        @timeit to "print iterations" begin
            println("iter = ", n, ", Δh_max = ", Δh_max, " at ", max_index)
        end

        # Update the old to the new solution
        @timeit to "update state" begin
            state.h .= state.hⁿ⁺¹
            state.u .= state.uⁿ⁺¹
            state.v .= state.vⁿ⁺¹
            state.w .= state.wⁿ⁺¹
            time += Δt
        end

    end

end

if isBenchmark
    # todo: how to get the physical array from the full array??
    h_max = find_maximum(state.hⁿ⁺¹)
    println(h_max)
else
    @timeit to "plot result" begin
        plot_model(grid, state)
        p = plot(bulge)
        display(p)
    end
end

end

print_timer(to)

end

function check_convergence!(state::State, solver_data::Solver_data)
    solver_data.Δh .= abs.(state.hⁿ⁺¹ - state.h)
    Δh_max, max_index = find_maximum(solver_data.Δh)
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
