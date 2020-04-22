# postprocessing.jl

using Plots

include("grids.jl")
include("model.jl")

# TODO: with multiple dispatch!

function prepare_plots(model_type, state)

    plotvars = []
    #model_type = model.model_type

    if model_type == "2DH"
        plotvars = Dict("h"=>state.h, "u"=>state.u, "v"=>state.v)
    elseif model_type == "2DV"
        plotvars = Dict("h"=>state.h, "u"=>state.u, "w"=>state.w)
    elseif model_type == "3D"
        plotvars = Dict("h"=>state.h, "u"=>state.u, "v"=>state.v, "w"=>state.w)
    end

    return plotvars
end

# TODO: with multiple dispatch!

function set_plot_properties(model_type)

    #model_type = model.model_type
    if model_type == "2DH"
        plot_title = ["water table h [m]" "horizontal velocity u [m/s]" "horizontal velocity v [m/s]" "vertical velocity w [m/s]"]
        plot_layout = (3,1)
    elseif model_type == "2DV"
        plot_title = ["water table h [m]" "horizontal velocity u [m/s]" "horizontal velocity v [m/s]" "vertical velocity w [m/s]"]
        plot_layout = (3,1)
    elseif model_type == "3D"
        plot_title = ["water table h [m]" "horizontal velocity u [m/s]" "horizontal velocity v [m/s]" "vertical velocity w [m/s]"]
        plot_layout = (4,1)
    end
    return plot_title, plot_layout
end

function plot_model(grid::Grid2dh, state::State2dh)
    # Plotting of 2DH model

    x = grid.x
    y = grid.y

    p1 = contour(x, y, state.h, fill = true)
    p2 = contour(x, y, state.u, fill = true)
    p3 = contour(x, y, state.v, fill = true)

    # For quvier plot of velocity
    velfac = 100000.0
    xvec = vec(grid.x)
    yvec = vec(grid.y)
    uvec = velfac.*vec(state.u)
    vvec = velfac.*vec(state.v)
    p4 = quiver(xvec,yvec, quiver=(uvec, vvec))
    display(plot(p1, p2, p3, p4, title=["water table h [m]" "horizontal velocity u [m/s]" "horizontal velocity v [m/s]" "velocity vectors [m/s]"], layout = @layout([a; b; c; d])))

end

function plot_model(state::State2dv)
    # Plotting of 2DV model

    x = grid.x
    z = grid.z

    p1 = contour(x, z, state.h, fill = true)
    p2 = contour(x, z, state.u, fill = true)
    p3 = contour(x, z, state.w, fill = true)
    display(plot(p1, p2, p3, title=["water table h [m]" "horizontal velocity u [m/s]" "vertical velocity w [m/s]"], layout = @layout([a; b; c])))

end

function plot_model(state::State3d)
    # Plotting of 3D model

    p1 = heatmap(state.h')
    p2 = heatmap(state.u')
    p3 = heatmap(state.v')
    p4 = heatmap(state.w')

    # Display the plots
    display(plot!(p1,p2,p3,p4, title=["water table h [m]" "horizontal velocity u [m/s]" "horizontal velocity v [m/s]" "vertical velocity w [m/s]"], layout=(4,1)))
end
