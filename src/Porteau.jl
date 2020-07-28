module Porteau

# Write your package code here.
include("athesis/operators.jl")
include("athesis/kernels_pressure_equation.jl")
include("athesis/kernels_darcy_equation.jl")

include("athesis/kernels.jl")
include("athesis/grids.jl")
include("athesis/model.jl")

include("athesis/gridloop.jl")
include("athesis/equations.jl")
include("athesis/model_input.jl")
include("athesis/external_forcing.jl")

include("athesis/postprocessing.jl")

include("athesis/initialize.jl")
include("athesis/reduce.jl")
include("athesis/reduce_max.jl")

export model_initialize
export set_sources!
export set_recharge!
export pressure_equation!
export darcy_equation!
export plot_model

end
