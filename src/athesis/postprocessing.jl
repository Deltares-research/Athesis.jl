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

    ilayer = grid.nz
    i_plot = Int(ceil(grid.nx/2))
    j_plot = Int(ceil(grid.ny/2))

    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    x = grid.x
    y = grid.y
    z = grid.z

    #p1 = heatmap(x,y,state.h[:,:,ilayer]', fill = true, c=cgrad([:white,:green,:blue]))
    p1 = plot(x[1:nx+2],y[1:ny+2],state.hⁿ⁺¹[1:nx+2,1:ny+2,ilayer], seriestype=:wireframe, camera=(30,60))
    #p2 = heatmap(state.u[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))
    #p3 = heatmap(state.v[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))
    #p4 = heatmap(state.w[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))
    p2 = plot(x[1:nx+2],state.hⁿ⁺¹[1:nx+2,j_plot,ilayer])
    #p3 = plot(x[:],state.uⁿ⁺¹[1:nx,Int(grid.ny/2),ilayer])
    p3 = plot(z[1:nz+2],state.wⁿ⁺¹[i_plot,j_plot,1:nz+2])
    # Display the plots
    display(plot(p1,p2,p3, layout = (3,1), size=(800,1200)))
    #display(plot!(p1,p2,p3,p4, title=["water table h [m]" "horizontal velocity u [m/s]" "horizontal velocity v [m/s]" "vertical velocity w [m/s]"], layout=(4,1)))
    #savefig("plot.png")
end
