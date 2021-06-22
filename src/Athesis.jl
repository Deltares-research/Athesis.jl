module Athesis

include("athesis/utils.jl")
include("athesis/initialize.jl")
include("athesis/run.jl")
include("athesis/equations.jl")
include("athesis/external_forcing.jl")
include("athesis/postprocessing.jl")
include("athesis/model.jl")
include("athesis/model_input.jl")
include("athesis/grids.jl")
include("athesis/boundary_conditions.jl")
include("athesis/gridloop.jl")
include("athesis/kernels_pressure_equation.jl")
include("athesis/kernels_darcy_equation.jl")
include("athesis/kernels.jl")
include("athesis/operators.jl")
include("athesis/fields.jl")
include("athesis/bmi.jl")


export getDefaultInput,
    initSimulation,
    runSimulation!,
    doExplicitIter!,
    plotSimulation,

    # Data structures
    Grid,
    Model,
    ModelInput,
    Parameters,
    State,
    TimeData,
    SolverData,
    Simulation,

    # Macros
    @synctimeit,
    @withCUDA

end
