# postprocessing.jl

using Plots

include("grids.jl")
include("model.jl")

# function prepare_plots(state)
#     plotvars = Dict("h"=>state.h, "u"=>state.u, "v"=>state.v, "w"=>state.w)
# end
#
# # TODO: with multiple dispatch!
#
# function set_plot_properties()
#     plot_title = ["water table h [m]" "horizontal velocity u [m/s]" "horizontal velocity v [m/s]" "vertical velocity w [m/s]"]
#     plot_layout = (4,1)
#     return plot_title, plot_layout
# end


function plot_model(grid, state)
    # Plotting of 3D model

    x = grid.x
    #y = grid.y
    z = grid.z

    nn = Int16(grid.nx/2)

    p1 = heatmap(x,z,state.h[nn,:,:]', fill = true, c=cgrad([:white,:green,:blue]))
    #p2 = heatmap(state.u[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))
    #p3 = heatmap(state.v[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))
    #p4 = heatmap(state.w[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))

    # Display the plots
    display(plot(p1))
    #display(plot!(p1,p2,p3,p4, title=["water table h [m]" "horizontal velocity u [m/s]" "horizontal velocity v [m/s]" "vertical velocity w [m/s]"], layout=(4,1)))
end
