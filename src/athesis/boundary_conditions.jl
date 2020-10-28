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
    value = bc.bc_pressure[1]
    for k = 2:nz+1
        for j = 2:ny+1
            h[1,j,k] = value
            hⁿ⁺¹[1,j,k] = value
        end
    end

    # East boundary
    value = bc.bc_pressure[2]
    for k = 2:nz+1
        for j = 2:ny+1
            h[nx+2,j,k] = value
            hⁿ⁺¹[nx+2,j,k] = value
        end
    end

    # Now the four Neumann boundary conditions (zero-flux)
    # at the south and north boundaries
    # and at the bottom and top boundaries

    # South boundary
    gradient = 0.0
    for k = 2:nz+1
        for i = 2:nx+1
            h[i,1,k] = h[i,2,k] - gradient*Δy
            hⁿ⁺¹[i,1,k] = h[i,1,k]
        end
    end

    # North boundary
    gradient = 0.0
    for k = 2:nz+1
        for i = 2:nx+1
            h[i,ny+2,k] = h[i,ny+1,k] + gradient*Δy
            hⁿ⁺¹[i,ny+2,k] = h[i,ny+2,k]
        end
    end

    # Bottom boundary
    gradient = 0.0
    for j = 2:ny+1
        for i = 2:nx+1
            h[i,j,1] = h[i,j,2] - gradient*Δz
            hⁿ⁺¹[i,j,1] = h[i,j,1]
        end
    end

    # Top boundary
    gradient = 0.0
    for j = 2:ny+1
        for i = 2:nx+1
            h[i,j,nz+2] = h[i,j,nz+1] + gradient*Δz
            hⁿ⁺¹[i,j,nz+2] = h[i,j,nz+2]
        end
    end

end
