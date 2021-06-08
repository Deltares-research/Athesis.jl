# '''
# Athesis.jl
# A Thramework for Efficient Simulations in Software
# A flexible and portable toolbox for scientific computing on CPUs and GPUs
# '''

using Test
using CuArrays
using CUDAdrv, CUDAnative

include("reduce.jl")

"""
Get model input:                                    
grid, model type, physcial and numerical parameters 
"""
function model_input(GPU)
    model_data = []
    grid_data  = []
    time_data  = []
    parameters = []
    err        = ' '

    println("Reading model input...")

    # Grid boundaries or south-west (lower-left) corner
    x0 = 0.0
    y0 = 0.0

    # grid dimensions
    nx = 101
    ny = 101

    # grid step size
    dx = 10.0
    dy = dx

    # Model components
    advection = false
    diffusion = true
    friction  = false

    # permeability?
    # 1.0e-06 m/s for some kind of sand
    # isotropic
    K = 1.0e-06

    # Courant number
    C = 0.98

    # time step size
    maxdx = max(dx, dy)
    dx2   = maxdx * maxdx
    dt    = C * 0.5 * dx2 / K

    # Maximum number of time steps
    maxsteps = 10000

    # Initialize hydraulic head
    h0 = 15.0 * ones(Float64, nx, ny)

    # boundary conditions
    bcleft = 10.0 * ones(Float64, nx)
    bcright = 20.0 * ones(Float64, nx)

    # solver parameters
    conv_limit = 1.0e-06

    # time integtation
    time_integration = "forward_euler"

    # Hardware/computing information
    # GPU = false

    grid_data  = Dict("x0" => x0, "y0" => y0, "nx" => nx, "ny" => ny, "dx" => dx, "dy" => dy)
    model_data = Dict("advection" => advection, "diffusion" => diffusion, "friction" => friction, "h0" => h0, "bcleft" => bcleft, "bcright" => bcright)
    time_data  = Dict("dt" => dt, "maxsteps" => maxsteps, "Courant" => C, "time_integration" => time_integration)
    parameters = Dict("K" => K, "conv_limit" => conv_limit, "GPU" => GPU)

    return model_data, grid_data, time_data, parameters, err
end


"""Setup the grid: structured / unstructured"""
function grid_setup(grid_data)

    println("Setting up the grid ...")
    grid = []
    err  = ' '

    x0 = grid_data["x0"]
    y0 = grid_data["y0"]
    nx = grid_data["nx"]
    ny = grid_data["ny"]
    dx = grid_data["dx"]
    dy = grid_data["dy"]

    xc = zeros(Float64, nx, ny)
    yc = zeros(Float64, nx, ny)
    for j = 1:ny
        for i = 1:nx
            xc[i,j] = x0 + (i - 0.5) * dx
            yc[i,j] = y0 + (j - 0.5) * dy
        end
    end

    grid = Dict("nx" => nx, "ny" => ny, "dx" => dx, "dy" => dy, "xc" => xc, "yc" => yc)

    return grid, err
end


"""
Setup the model:                                  
discretization in space and time, initializations 
"""
function model_setup(model_data, time_data, grid, parameters)

    println("Setting up the model ...")

    model = []
    pars  = []
    err   = ' '

    advection = model_data["advection"]
    diffusion = model_data["diffusion"]
    friction  = model_data["friction"]
    h0        = model_data["h0"]
    bcleft    = model_data["bcleft"]
    bcright   = model_data["bcright"]

    bcs       = [bcleft, bcright]

    dt               = time_data["dt"]
    maxsteps         = time_data["maxsteps"]
    C                = time_data["Courant"]
    # time_integration = time_data["time_integration"]

    conv_limit = parameters["conv_limit"]
    GPU        = parameters["GPU"]

    nx = grid["nx"]
    ny = grid["ny"]
    # dx = grid["dx"]
    # dy = grid["dy"]

    # Set initial h
    h = h0

    # Special GPU initializations and copying
    if GPU
        # copy to GPU
        h = CuArray(h)
        hnew = copy(h)
        hchange = 2.0 * conv_limit * CuArrays.ones(eltype(h), nx * ny)
        state = [h, hnew, hchange]

        # GPU stuff
        ths = (16, 16)
        bls = (Int(ceil(nx / ths[1])), Int(ceil(ny / ths[2])))
        pars  = [ths, bls]
    else
        hnew = copy(h)
        hchange = 2.0 * conv_limit * ones(eltype(h), nx * ny)
    end

    state = [h, hnew, hchange]
    model = Dict("state" => state, "bcs" => bcs)

    # For now the simple diffusion model, with forward euler (explicit) time integration
    # TO DO

    return model, pars, err
end

"""Kernel for change in 2D, result in u"""
function calculateChange_gpu!(u::AbstractArray, v::AbstractArray, w::AbstractArray)
    ix = blockDim().x * (blockIdx().x - 1) + threadIdx().x;
    iy = blockDim().y * (blockIdx().y - 1) + threadIdx().y;
    idx = size(u, 1) * (iy - 1) + ix
    if ix < size(u, 1) + 1 && iy < size(u, 2) + 1
        @inbounds w[idx] = abs(v[ix,iy] - u[ix,iy])
    end

    return nothing
end

"""Kernel for diffusion in 2D"""
function calculateDiffusion_gpu!(unew::AbstractArray, u::AbstractArray, dx, dy, dt, K)
    ix = blockDim().x * (blockIdx().x - 1) + threadIdx().x;
    iy = blockDim().y * (blockIdx().y - 1) + threadIdx().y;

    if ix > 1 && iy > 1 && ix < size(unew, 1) && iy < size(unew, 2)
        @inbounds F = K * ( (u[ix + 1,iy] + u[ix - 1,iy] - 2.0 * u[ix,iy]) / (dx * dx) + (u[ix,iy + 1] + u[ix,iy - 1] - 2.0 * u[ix,iy]) / (dy * dy) )
        @inbounds unew[ix,iy] = u[ix,iy] + dt * F
    end

    return nothing
end

function update_boundary_conditions!(model)
    bcs = []
    err = ' '

    # Boundary conditions
    bcleft    = model["bcs"][1]
    bcright   = model["bcs"][2]

    # Unpack present state
    h         = model["state"][1]
    hnew      = model["state"][2]
    hchange   = model["state"][3]

    # set boundaries
    h[:,1] = bcleft
    h[:,end] = bcright

    state = [h,hnew, hchange]
    model["state"] = state

    return model, err
end


"""Time loop: integrate the solution in time"""
function timeloop!(model, pars, grid, time_data, parameters)

    println("Starting time loop ...")
    # solution = []
    err      = ' '

    dt       = time_data["dt"]
    maxsteps = time_data["maxsteps"]
    # state    = model["state"]
    GPU      = parameters["GPU"]

    # niter = 0
    conv = false

    for n = 1:maxsteps
        if GPU
            model, conv, err = single_timestep_GPU!(n, dt, model, pars, grid, parameters, conv)
        else
            model, conv, err = single_timestep!(n, dt, model, grid, parameters, conv)
        end

        if conv
            break
        end
    end

    # println("convergence in ", niter, " iterations (", maxchange, ")")
    # display(hnew)
    display(model["state"][2])

    return model, err
end

"""Perform a single time step on the GPU"""
function single_timestep_GPU!(niter, dt, model, pars, grid, parameters, conv)
    err = ' '
    maxchange = 0.0

    # set boundaries for present time steps
    model, err = update_boundary_conditions!(model)

    # Unpack present state
    h       = model["state"][1]
    hnew    = model["state"][2]
    hchange = model["state"][3]

    conv_limit = parameters["conv_limit"]
    K          = parameters["K"]

    # GPU parameters
    ths = pars[1]
    bls = pars[2]

    # Unpack grid
    dx = grid["dx"]
    dy = grid["dy"]
    # nx = grid["nx"]
    # ny = grid["ny"]

    # niter = niter + 1

    # run CUDA kernels
    @cuda threads = ths blocks = bls calculateDiffusion_gpu!(hnew, h, dx, dy, dt, K)
    @cuda threads = ths blocks = bls calculateChange_gpu!(h, hnew, hchange)

    # reduce to find maximum change in h
    gpu_reduce(max, hchange, hchange)
    maxchange = hchange[1]

    # set no-flow boundaries
    hnew[1,:] = hnew[2,:]
    hnew[end,:] = hnew[end - 1,:]

    # swap
    htemp = hnew
    hnew  = h
    h     = htemp

    state = [h, hnew, hchange]

    if maxchange < conv_limit
        conv = true
    end
    # maxchange = 0.0

    if (mod(niter, 10) == 0)
        println("iteration: ", niter, ". Maximum change: ", maxchange)
    end
    model["state"] = state

    return model, conv, err
end

"""Perform a single time step"""
function single_timestep!(niter, dt, model, grid, parameters, conv)
    err = ' '
    maxchange = 0.0

    # set boundaries for present time steps
    model, err = update_boundary_conditions!(model)

    # Unpack present state
    h       = model["state"][1]
    hnew    = model["state"][2]
    hchange = model["state"][3]

    conv_limit = parameters["conv_limit"]
    K = parameters["K"]

    dx = grid["dx"]
    dy = grid["dy"]
    nx = grid["nx"]
    ny = grid["ny"]
    # niter = niter + 1

    # Compute the spatial operator
    for j = 2:ny - 1
        for i = 2:nx - 1
            @inbounds F = K * ( (h[i + 1,j] + h[i - 1,j] - 2.0 * h[i,j]) / (dx * dx) + (h[i,j + 1] + h[i,j - 1] - 2.0 * h[i,j]) / (dy * dy) )
            @inbounds hnew[i,j] = h[i,j] + dt * F

            change = abs(h[i,j] - hnew[i,j])
            # replace by max(..)?
            if change > maxchange
                maxchange = change
            end
        end
    end

    # set no-flow boundaries
    hnew[1,:] = hnew[2,:]
    hnew[end,:] = hnew[end - 1,:]

    # swap
    htemp = hnew
    hnew  = h
    h     = htemp

    model["state"] = [h, hnew, hchange]

    if maxchange < conv_limit
        conv = true
    end
    # maxchange = 0.0
    if (mod(niter, 10) == 0)
        println("iteration: ", niter, ". Maximum change: ", maxchange)
    end

    # solution = state

    return model, conv, err
end

"""Finalize the program"""
function finalization!(model, grid)

    println("Finalizing model ...")
    err = ' '

    model = nothing
    grid  = nothing
    return err
end

"""Main program"""
function athesis(; GPU::Bool=false)

    println("Running Athesis ...")

    # Model input
    model_data, grid_data, time_data, parameters, err = model_input(GPU)

    # Grid setup
    grid, err = grid_setup(grid_data)

    # Model setup
    model, pars, err = model_setup(model_data, time_data, grid, parameters)

    # Start of the time loop
    solution, err = timeloop!(model, pars, grid, time_data, parameters)

    # Finalization
    err = finalization!(model, grid)

    return nothing

end

function time_athesis(; GPU::Bool=false)

    @time athesis(GPU=GPU)

end

time_athesis(GPU=true)
