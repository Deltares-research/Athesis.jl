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