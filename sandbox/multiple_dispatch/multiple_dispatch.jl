# multiple_dispatch test

ENV["MPLBACKEND"]="tkagg"
using PyPlot
pygui(true)

mutable struct Grid
    nx::Int64
    dx::Float64
end

mutable struct Params
    K::Float64
    dt::Float64
    nsteps::Int64
    Tend::Float64
    Courant::Float64
end

function preprocess(method)
    if method == "advection_explicit"
        advection       = true
        advection_theta = false
        diffusion       = 0
        diffusion_theta = false
    elseif method == "advection_implicit"
        advection       = true
        advection_theta = 1.0
        diffusion       = 0
        diffusion_theta = false
    elseif method == "diffusion_explicit"
        advection       = 0
        advection_theta = false
        diffusion       = true
        diffusion_theta = false
    elseif method == "diffusion_implicit"
        advection_theta = false
        advection       = 0
        diffusion_theta = 1.0
        diffusion       = true
    elseif method == "advection_explicit_diffusion_explicit"
        advection       = true
        advection_theta = false
        diffusion       = true
        diffusion_theta = false
    elseif method == "advection_implicit_diffusion_explicit"
        advection       = true
        advection_theta = 1.0
        diffusion       = true
        diffusion_theta = false
    elseif method == "advection_explicit_diffusion_implicit"
        advection       = true
        advection_theta = false
        diffusion       = true
        diffusion_theta = 1.0
    elseif method == "advection_implicit_diffusion_implicit"
        advection       = true
        advection_theta = 1.0
        diffusion       = true
        diffusion_theta = 1.0
    end

    #model_components = [advection, advection_theta, diffusion, diffusion_theta]
    #return model_components
    return advection, advection_theta, diffusion, diffusion_theta

end

function initialize!(model_components)

    L      = 1000
    nx     = 1000
    dx     = L/nx
    K      = 0.1
    nsteps = 1000
    dt0    = 5.0
    Tend   = 100000.0

    # Courant number
    Courant = 0.98
    maxdx  = dx #max(dx,dy)

    # Unpack model_components
    advection       = model_components[1]
    advection_theta = model_components[2]
    diffusion       = model_components[3]
    diffusion_theta = model_components[4]

    # time step size based on diffusion
    if diffusion==true && diffusion_theta == false
        dx2   = maxdx*maxdx
        dt0    = min(dt0, Courant*0.5*dx2/K)
    end
    dt = dt0

    params = Params(K, dt, nsteps, Tend, Courant)

    h    = Array{Float64,1}(undef, nx)
    diag = 1.0*ones(Float64, nx)
    rhs  = Array{Float64,1}(undef, nx)

    # Save state
    state = [h, diag, rhs]

    # Set initial solution
    for i = 1:nx
        if nx*0.4 < i < nx*0.6
            h[i] = 0.4
        else
            h[i] = 0.2
        end
    end
    initial_condition = copy(h)

    grid   = Grid(nx, dx)

    return grid, state, initial_condition, params
end

function compdt!(params, model_components, grid, time, initial_condition)
    # time step size based on advection

    # Unpack model_components
    advection       = model_components[1]
    advection_theta = model_components[2]
    diffusion       = model_components[3]
    diffusion_theta = model_components[4]

    # Local variable dt
    #dt = 0.0

    if advection==true && advection_theta==false
        hmax = maximum(abs.(initial_condition))
        params.dt    = min(params.dt, params.Courant*grid.dx/hmax)

        # Don't exceed end time
        params.dt    = min(params.dt, params.Tend-time)
    end
    #println(dt)
    return params.dt
end

# Update function for explicit advection
function update!(state, grid, params, adv::Bool, diff::Int64, adv_theta::Bool, diff_theta::Bool)
    # Implementation corresponding to method:
    # explicit advection
    #println("Update state using explicit advection")
    state = advection_explicit!(state, grid, params)
    return state
end

# Update function for implicit advection
function update!(state, grid, params, adv::Bool, diff::Int64, adv_theta::Float64, diff_theta::Bool)
    # Implementation corresponding to method:
    # implicit advection
    #println("Update state using implicit advection")
    state = advection_implicit!(state, grid, params, adv_theta)
    return state
end

# Update function for explicit diffusion
function update!(state, grid, params, adv::Int64, diff::Bool, adv_theta::Bool, diff_theta::Bool)
    # Implementation corresponding to method:
    # explicit diffusion
    #println("Update state using explicit diffusion")
    state = diffusion_explicit!(state, grid, params)
    return state
end

# Update function for implicit diffusion
function update!(state, grid, params, adv::Int64, diff::Bool, adv_theta::Bool, diff_theta::Float64)
    # Implementation corresponding to method:
    # implicit diffusion
    #println("Update state using implicit diffusion")
    state = diffusion_implicit!(state, grid, params, diff_theta)
    return state
end

# Update function for explicit advection and explicit diffusion
function update!(state, grid, params, adv::Bool, diff::Bool, adv_theta::Bool, diff_theta::Bool)
    # Implementation corresponding to method:
    # explicit advection and explicit diffusion
    #println("Update state using explicit advection and explicit diffusion")
    state = advection_explicit!(state, grid, params)
    state = diffusion_explicit!(state, grid, params)
    return state
end

# Update function for implicit advection and explicit diffusion
function update!(state, grid, params, adv::Bool, diff::Bool, adv_theta::Float64, diff_theta::Bool)
    # Implementation corresponding to method:
    # implicit advection and explicit diffusion
    #println("Update state using implicit advection and explicit diffusion")
    state = advection_implicit!(state, grid, params, adv_theta)
    state = diffusion_explicit!(state, grid, params)
    return state
end

# Update function for explicit advection and implicit diffusion
function update!(state, grid, params, adv::Bool, diff::Bool, adv_theta::Bool, diff_theta::Float64)
    # Implementation corresponding to method:
    # explicit advection and implicit diffusion
    #println("Update state using explicit advection and implicit diffusion")
    state = advection_explicit!(state, grid, params)
    state = diffusion_implicit!(state, grid, params, diff_theta)
    return state
end

# Update function for implicit advection and implicit diffusion
function update!(state, grid, params, adv::Bool, diff::Bool, adv_theta::Float64, diff_theta::Float64)
    # Implementation corresponding to method:
    # implicit advection and implicit diffusion
    #println("Update state using implicit advection and implicit diffusion")
    state = advection_implicit!(state, grid, params, adv_theta)
    state = diffusion_implicit!(state, grid, params, diff_theta)
    return state
end

# Explicit advection computation
function advection_explicit!(state, grid, params)

    nx  = grid.nx
    dx  = grid.dx

    dt  = params.dt

    h    = state[1]
    diag = state[2]
    rhs  = state[3]

    for i=2:nx-1
        uxm = (h[i]-h[i-1])/dx
        uxp = (h[i+1]-h[i])/dx
        am  = min(h[i],0.0)
        ap  = max(h[i],0.0)
        adv = ap*uxm + am*uxp
        @inbounds rhs[i] +=  -dt * adv
    end

    state = [h, diag, rhs]

    return state
end

# Implicit advection computation
function advection_implicit!(state, grid, params, adv_theta)
    nx  = grid.nx
    dx  = grid.dx

    dt  = params.dt

    h    = state[1]
    diag = state[2]
    rhs  = state[3]

    for i=2:nx-1
        uxm = -h[i-1]/dx
        uxp = h[i+1]/dx
        am  = min(h[i],0.0)
        ap  = max(h[i],0.0)
        adv_expl = ap*uxm + am*uxp
        adv_impl = ap/dx - am/dx
        @inbounds rhs[i]  += -dt * adv_expl
        @inbounds diag[i] +=  dt * adv_impl
    end

    state = [h, diag, rhs]

    return state
end

# Explicit advection computation
function diffusion_explicit!(state, grid, params)

    nx  = grid.nx
    dx  = grid.dx
    dx2 = dx*dx

    K   = params.K
    dt  = params.dt
    dtK = dt*K

    h    = state[1]
    diag = state[2]
    rhs  = state[3]

    #hnew = copy(h)

    for i=2:nx-1
        @inbounds rhs[i] += dtK * (h[i+1] - 2.0*h[i] + h[i-1]) / dx2
        #@inbounds hnew[i,j] = h[i,j] + dt*K * ( (h[i+1,j] + h[i-1,j] - 2.0*h[i,j])/(dx*dx) + (h[i,j+1]+h[i,j-1] - 2.0*h[i,j])/(dy*dy) )
    end

    state = [h, diag, rhs]

    return state

end

# Implicit diffusion computation
function diffusion_implicit!(state, grid, params, diff_theta)
    nx  = grid.nx
    dx  = grid.dx
    dx2 = dx*dx

    K   = params.K
    dt  = params.dt
    dtK = dt*K

    h    = state[1]
    diag = state[2]
    rhs  = state[3]

    #hnew = copy(h)

    for i=2:nx-1
        @inbounds diag[i] += dtK * 2.0 / dx2
        @inbounds rhs[i] += dtK * (h[i+1] + h[i-1]) / dx2
        #@inbounds hnew[i,j] = h[i,j] + dt*K * ( (h[i+1,j] + h[i-1,j] - 2.0*h[i,j])/(dx*dx) + (h[i,j+1]+h[i,j-1] - 2.0*h[i,j])/(dy*dy) )
    end

    state = [h, diag, rhs]

    return state
end

# Functio for solving a diagonal system of equations
function solve!(state, grid)
    h    = state[1]
    diag = state[2]
    rhs  = state[3]
    nx   = grid.nx

    for i = 1:nx
        h[i] = rhs[i]/diag[i]
        diag[i] = 1.0
    end

    return state
end

function plot_solution(n, time, h)
    clf()
    display(plot(h, linestyle="-",marker="o"))
    title("Time = $time s. (timestep $n)")
    pause(0.1)
end

# Test multiple dispatch for update function
# based on the contents of "parameters"
# For now for a 1D model
# depending on the argument, 8 methods can be invoked
# (through multiple dispatch of UPDATE function):
#   - explicit advection
#   - explicit diffusion
#   - (diagonally-)implicit advection
#   - (diagonally-)implicit diffusion
#   - explicit advection, explicit diffusion
#   - (diagonally-)implicit advection, explicit diffusion
#   - (diagonally-)implicit advection, explicit diffusion
#   - (diagonally-)implicit advection, (diagonally-)implicit diffusion

function test_multiple_dispatch(;method="advection_explicit")

    # Preprocess function arguments to determine model components
    # needed for multiple dispatch of UPDATE function
    advection, advection_theta, diffusion, diffusion_theta = preprocess(method)

    # Pack model components
    model_components = [advection, advection_theta, diffusion, diffusion_theta]

    # Initialize model
    grid, state, initial_condition, params = initialize!(model_components)

    # Number of grid points in x-direction
    nx = grid.nx

    # Unpack state
    h    = state[1]
    diag = state[2]
    rhs  = state[3]

    time = 0.0
    for n = 1:params.nsteps

        # Determine time step (possible time step restriction from explicit model components)
        params.dt     = compdt!(params, model_components, grid, time, initial_condition)

        # Update time
        time  += params.dt

        # Update state
        rhs    = copy(h)
        state  = [h, diag, rhs]
        state  = update!(state, grid, params, advection, diffusion, advection_theta, diffusion_theta)

        # Solve system of equations (diagonal for now)
        state  = solve!(state, grid)

        # Plot solution
        if mod(n,50) == 0
            plot_solution(n, time, h)
        end

        # Exit at end time
        if time == params.Tend
            break
        end
    end

    close()
    println("\n")
end

#@time test_multiple_dispatch(;method="advection_explicit")
#@time test_multiple_dispatch(;method="advection_implicit")
#@time test_multiple_dispatch(;method="diffusion_explicit")
#@time test_multiple_dispatch(;method="diffusion_implicit")
#@time test_multiple_dispatch(;method="advection_explicit_diffusion_explicit")
#@time test_multiple_dispatch(;method="advection_implicit_diffusion_explicit")
#@time test_multiple_dispatch(;method="advection_explicit_diffusion_implicit")
@time test_multiple_dispatch(;method="advection_implicit_diffusion_implicit")
