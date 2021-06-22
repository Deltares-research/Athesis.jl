# Basic Model Interface (BMI) implementation based on
# https://github.com/Deltares/BasicModelInterface.jl

using BasicModelInterface
const BMI = BasicModelInterface


"""
    BMI.initialize(::Type{<:Athesis.Simulation}, config_file)

Initialize the simulation.
Will return a Simulation that is ready to run.
"""
function BMI.initialize(::Type{<:Simulation}, config_file)
    myFloat = Float64
    input = getDefaultInput(myFloat)
    useGPU = false
    simulation = initSimulation(input, useGPU, myFloat)
    return simulation
end

"""
    BMI.update(simulation::Simulation)

Advance model state by one time step.

Perform all tasks that take place within one pass through the model's
time loop.
"""
function BMI.update(simulation::Simulation)
    # We do not know the number of timesteps so set it to 1
    n = 1
    doTimestep!(n, simulation)
end

"""
    BMI.update_until(simulation::Simulation)

Advance model state until the given time.

The given `time` must be a model time later than the current model time.
"""
function BMI.update_until(simulation::Simulation, time)
    # We do not know the number of timesteps so set it to 1
    n = 1
    while simulation.timeData.time < time
        doTimestep!(n, simulation)
    end
end

"""
    BMI.get_start_time(simulation::Simulation)

Start time of the simulation.
"""
function BMI.get_start_time(simulation::Simulation)
    return 0.0
end

"""
    BMI.get_end_time(simulation::Simulation)

End time of the simulation.
"""
function BMI.get_end_time(simulation::Simulation)
    return simulation.timeData.tend
end

"""
    BMI.get_current_time(simulation::Simulation)

Current time of the simulation.
"""
function BMI.get_current_time(simulation::Simulation)
    return simulation.timeData.time
end

"""
    BMI.get_time_step(simulation::Simulation)

Time step of the simulation.
"""
function BMI.get_time_step(simulation::Simulation)
    return simulation.timeData.Δt
end