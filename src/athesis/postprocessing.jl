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

    ilayer = 2

    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    x = grid.x
    y = grid.y
    #p1 = heatmap(x,y,state.h[:,:,ilayer]', fill = true, c=cgrad([:white,:green,:blue]))
    p1 = plot(x,y,state.hⁿ⁺¹[1:nx,1:ny,ilayer], seriestype=:wireframe, camera=(30,60))
    #p2 = heatmap(state.u[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))
    #p3 = heatmap(state.v[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))
    #p4 = heatmap(state.w[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))
    p2 = plot(x[:],state.hⁿ⁺¹[1:nx,Int(grid.ny/2),ilayer])
    # Display the plots
    display(plot(p1,p2, layout = (2,1), size=(800,800)))
    #display(plot!(p1,p2,p3,p4, title=["water table h [m]" "horizontal velocity u [m/s]" "horizontal velocity v [m/s]" "vertical velocity w [m/s]"], layout=(4,1)))
    #savefig("plot.png")
end
