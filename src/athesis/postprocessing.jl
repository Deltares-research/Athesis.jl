# postprocessing.jl

using Plots

function plotSimulation(simulation, to::TimerOutput)
    # Plotting of 3D model
    @synctimeit to "plot result" begin
        grid = simulation.grid
        state = simulation.state


        ilayer = grid.nz

        nx = grid.nx
        ny = grid.ny
        nz = grid.nz

        x = grid.x
        y = grid.y
        z = grid.z

        # p1 = heatmap(x,y,state.h[:,:,ilayer]', fill = true, c=cgrad([:white,:green,:blue]))
        p1 = plot(x[0:nx + 1], y[0:ny + 1], state.hⁿ⁺¹[0:nx + 1,0:ny + 1,ilayer], seriestype=:wireframe, camera=(30, 60))
        # p2 = heatmap(state.u[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))
        # p3 = heatmap(state.v[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))
        # p4 = heatmap(state.w[:,nn,:]', fill = true, c=cgrad([:white,:green,:blue]))
        p2 = plot(x[0:nx + 1], state.hⁿ⁺¹[0:nx + 1,Int(ceil(grid.ny / 2)),ilayer])
        # p3 = plot(x[:],state.uⁿ⁺¹[1:nx,Int(grid.ny/2),ilayer])
        p3 = plot(z[0:nz + 1], state.wⁿ⁺¹[Int(ceil(grid.nx / 2)),Int(ceil(grid.ny / 2)),0:nz + 1])
        # Display the plots
        display(plot(p1, p2, p3, layout=(3, 1), size=(800, 1200)))
        # display(plot!(p1,p2,p3,p4, title=["water table h [m]" "horizontal velocity u [m/s]" "horizontal velocity v [m/s]" "vertical velocity w [m/s]"], layout=(4,1)))
        # savefig("plot.png")
        
    end # end timer
end
