# Boundary conditions

function set_boundaries!(grid, state, bc)
    # For now we simply set 2 Dirichlet boundary conditions on the west and east boundaries
    # and we set 4 Neumann boundary conditions on the south, north, bottom and top boundaries
    # All boundaries are pressure boundaries for now

    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    Δx = grid.Δx
    Δy = grid.Δy
    Δz = grid.Δz

    h = state.h
    hⁿ⁺¹ = state.hⁿ⁺¹

    # West  boundary
    value = bc.west_pressure
    for k = 1:nz
        for j = 1:ny
            h[0,j,k] = value
            hⁿ⁺¹[0,j,k] = value
        end
    end

    # East boundary
    value = bc.east_pressure
    for k = 1:nz
        for j = 1:ny
            h[nx+1,j,k] = value
            hⁿ⁺¹[nx+1,j,k] = value
        end
    end

    # Now the four Neumann boundary conditions (zero-flux)
    # at the south and north boundaries
    # and at the bottom and top boundaries

    # South boundary
    gradient = 0.0
    for k = 1:nz
        for i = 1:nx
            h[i,0,k] = h[i,1,k] - gradient*Δy
            hⁿ⁺¹[i,0,k] = h[i,0,k]
        end
    end

    # North boundary
    gradient = 0.0
    for k = 1:nz
        for i = 1:nx
            h[i,ny+1,k] = h[i,ny,k] + gradient*Δy
            hⁿ⁺¹[i,ny+1,k] = h[i,ny+1,k]
        end
    end

    # Bottom boundary
    gradient = 0.0
    for j = 1:ny
        for i = 1:nx
            h[i,j,0] = h[i,j,1] - gradient*Δz
            hⁿ⁺¹[i,j,0] = h[i,j,0]
        end
    end

    # Top boundary
    gradient = 0.0
    for j = 1:ny
        for i = 1:nx
            h[i,j,nz+1] = h[i,j,nz] + gradient*Δz
            hⁿ⁺¹[i,j,nz+1] = h[i,j,nz+1]
        end
    end

end
