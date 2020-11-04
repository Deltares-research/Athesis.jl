module Athesis

include("athesis/initialize.jl")
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
include("athesis/operators.jl")

export
    model_initialize,
    set_sources!,
    set_recharge!,
    set_boundaries!,
    pressure_equation!,
    darcy_equation!,
    plot_model,

    Grid,
    Model,
    Parameters,
    State,
    Time_data,
    Solver_data

end
