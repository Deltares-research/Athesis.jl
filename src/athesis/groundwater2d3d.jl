# 3D groundwater

using Plots
using CuArrays
using TimerOutputs

include("initialize.jl")
include("equations.jl")
include("external_forcing.jl")
include("postprocessing.jl")
include("model.jl")
include("boundary_conditions.jl")

# This is the present data storage:
# grid       = (nx, ny, nz, Δx, Δy, Δz, x, y, z)
# state      = (h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
# model      = externals, recharge
# parameters = (K, Δt, tend)
# externals  = source
# source     = (i_src, j_src, k_src, n_src, externals)


function groundwater3d()

    to = TimerOutput()

    @timeit to "run time" begin

    # Initialize the model
    println("Running 3D groundwater model:")

    @timeit to "initialization" begin
        grid, model, state, parameters, time_data, solver_data = model_initialize()

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
            #bc = model.boundary_conditions
            #new2old!(state,bc)
            state.h = copy(state.hⁿ⁺¹)
            state.u = copy(state.uⁿ⁺¹)
            state.v = copy(state.vⁿ⁺¹)
            state.w = copy(state.wⁿ⁺¹)
            time += Δt
        end

    end

end

@timeit to "plot result" begin
    plot_model(grid, state)
    p = plot(bulge)
    display(p)
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


function new2old!(state::State,bc::BoundaryConditions)
    h    = state.h
    u    = state.u
    v    = state.v
    w    = state.w
    hⁿ⁺¹ = state.hⁿ⁺¹
    uⁿ⁺¹ = state.uⁿ⁺¹
    vⁿ⁺¹ = state.vⁿ⁺¹
    wⁿ⁺¹ = state.wⁿ⁺¹

    bcpres = bc.bc_pressure

    copy_new2old!(bcpres, h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
end

function copy_new2old!(bcpres::Array, h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
    
    #println("CPU copy")
    for k = 0:size(h)[3]-1
        for j = 0:size(h)[2]-1
            for i = 0:size(h)[1]-1
                copy_kernel!(i, j, k,
                       h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
            end
        end
    end
end

function copy_new2old!(bcpres::CuArray, h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)

    #println("GPU copy")
    ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
    iy = (blockIdx().y-1)*blockDim().y  + threadIdx().y
    iz = (blockIdx().z-1)*blockDim().z  + threadIdx().z

    ths = (8,8,4)
    nbx = Int(ceil(size(h)[1]/ths[1]))
    nby = Int(ceil(size(h)[2]/ths[2]))
    nbz = Int(ceil(size(h)[3]/ths[3]))
    bls = (nbx,nby,nbz)

    @cuda threads=ths blocks=bls cuda_wrap_copy!(
                            h, u, v, w,
                            hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)

end

function cuda_wrap_copy!(h, u, v, w,
                         hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)

    ix = (blockIdx().x-1)*blockDim().x  + threadIdx().x
    iy = (blockIdx().y-1)*blockDim().y  + threadIdx().y
    iz = (blockIdx().z-1)*blockDim().z  + threadIdx().z
    if (0 <= ix <= nx+1 && 0 <= iy <= ny+1 && 0 <= iz <= nz+1)
        copy_kernel!(ix, iy, iz,
               h, u, v, w,
               hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹)
    end

    return nothing
end

function copy_kernel!(i,j,k,
                      h, u, v, w,
                      hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹ )

    h[i,j,k] = hⁿ⁺¹[i,j,k]
    u[i,j,k] = uⁿ⁺¹[i,j,k]
    v[i,j,k] = vⁿ⁺¹[i,j,k]
    w[i,j,k] = wⁿ⁺¹[i,j,k]
end

################################
# Execute the main function

@time groundwater3d()
################################
