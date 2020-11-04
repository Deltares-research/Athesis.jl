# test stencil operator

include("grids.jl")
#include("operators.jl")

abstract type Stencil end

# Internal cell type
struct InternalCell <:Stencil end

# 6 faces of the domain
struct WestBoundaryCell <:Stencil end
struct EastBoundaryCell <:Stencil end
struct SouthBoundaryCell <:Stencil end
struct NorthBoundaryCell <:Stencil end
struct BottomBoundaryCell <:Stencil end
struct TopBoundaryCell <:Stencil end

# 12 sides of the domain
struct SouthWestBoundaryCell <:Stencil end
struct SouthEastBoundaryCell <:Stencil end
struct NorthWestBoundaryCell <:Stencil end
struct NorthEastBoundaryCell <:Stencil end
struct BottomWestBoundaryCell <:Stencil end
struct BottomEastBoundaryCell <:Stencil end
struct TopWestBoundaryCell <:Stencil end
struct TopEastBoundaryCell <:Stencil end
struct BottomSouthBoundaryCell <:Stencil end
struct BottomNorthBoundaryCell <:Stencil end
struct TopSouthBoundaryCell <:Stencil end
struct TopNorthBoundaryCell <:Stencil end

# 8 corners of the domain
struct BottomSouthWestBoundaryCell <:Stencil end
struct BottomSouthEastBoundaryCell <:Stencil end
struct BottomNorthWestBoundaryCell <:Stencil end
struct BottomNorthEastBoundaryCell <:Stencil end
struct TopSouthWestBoundaryCell <:Stencil end
struct TopSouthEastBoundaryCell <:Stencil end
struct TopNorthWestBoundaryCell <:Stencil end
struct TopNorthEastBoundaryCell <:Stencil end


# Stencil for internal cell
function get_stencil!(cell, stencil_type::InternalCell, stencil, nb)
    (i, j, k) = cell
    nb = 6
    #stencil = Array{Tuple,1}(undef,6)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i+1, j  , k  )
    stencil[3] = (i  , j-1, k  )
    stencil[4] = (i  , j+1, k  )
    stencil[5] = (i  , j  , k-1)
    stencil[6] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for west boundary cell
function get_stencil!(cell, stencil_type::WestBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 5
    #stencil = Array{Tuple,1}(undef,5)
    stencil[1] = (i+1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j+1, k  )
    stencil[4] = (i  , j  , k-1)
    stencil[5] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for east boundary cell
function get_stencil!(cell, stencil_type::EastBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 5
    #stencil = Array{Tuple,1}(undef,5)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j+1, k  )
    stencil[4] = (i  , j  , k-1)
    stencil[5] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for south boundary cell
function get_stencil!(cell, stencil_type::SouthBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 5
    #stencil = Array{Tuple,1}(undef,5)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i+1, j  , k  )
    stencil[3] = (i  , j+1, k  )
    stencil[4] = (i  , j  , k-1)
    stencil[5] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for north boundary cell
function get_stencil!(cell, stencil_type::NorthBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 5
    #stencil = Array{Tuple,1}(undef,5)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i+1, j  , k  )
    stencil[3] = (i  , j-1, k  )
    stencil[4] = (i  , j  , k-1)
    stencil[5] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for bottom boundary cell
function get_stencil!(cell, stencil_type::BottomBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 5
    #stencil = Array{Tuple,1}(undef,5)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i+1, j  , k  )
    stencil[3] = (i  , j-1, k  )
    stencil[4] = (i  , j+1, k  )
    stencil[5] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for top boundary cell
function get_stencil!(cell, stencil_type::TopBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 5
    #stencil = Array{Tuple,1}(undef,5)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i+1, j  , k  )
    stencil[3] = (i  , j-1, k  )
    stencil[4] = (i  , j+1, k  )
    stencil[5] = (i  , j  , k-1)
    return stencil, nb
end

# Stencil for south-west boundary cell
function get_stencil!(cell, stencil_type::SouthWestBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i+1, j  , k  )
    stencil[2] = (i  , j+1, k  )
    stencil[3] = (i  , j  , k-1)
    stencil[4] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for south-west boundary cell
function get_stencil!(cell, stencil_type::SouthEastBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i  , j+1, k  )
    stencil[3] = (i  , j  , k-1)
    stencil[4] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for south-west boundary cell
function get_stencil!(cell, stencil_type::NorthWestBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i+1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j  , k-1)
    stencil[4] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for south-west boundary cell
function get_stencil!(cell, stencil_type::NorthEastBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j  , k-1)
    stencil[4] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for bottom-west boundary cell
function get_stencil!(cell, stencil_type::BottomWestBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i+1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j+1, k  )
    stencil[4] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for bottom-east boundary cell
function get_stencil!(cell, stencil_type::BottomEastBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j+1, k  )
    stencil[4] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for bottom-south boundary cell
function get_stencil!(cell, stencil_type::BottomSouthBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i+1, j  , k  )
    stencil[3] = (i  , j+1, k  )
    stencil[4] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for bottom-north boundary cell
function get_stencil!(cell, stencil_type::BottomNorthBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i+1, j  , k  )
    stencil[3] = (i  , j-1, k  )
    stencil[4] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for top-west boundary cell
function get_stencil!(cell, stencil_type::TopWestBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i+1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j+1, k  )
    stencil[4] = (i  , j  , k-1)
    return stencil, nb
end

# Stencil for top-east boundary cell
function get_stencil!(cell, stencil_type::TopEastBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j+1, k  )
    stencil[4] = (i  , j  , k-1)
    return stencil, nb
end

# Stencil for top-south boundary cell
function get_stencil!(cell, stencil_type::TopSouthBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i+1, j  , k  )
    stencil[3] = (i  , j+1, k  )
    stencil[4] = (i  , j  , k-1)
    return stencil, nb
end

# Stencil for top-north boundary cell
function get_stencil!(cell, stencil_type::TopNorthBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 4
    #stencil = Array{Tuple,1}(undef,4)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i+1, j  , k  )
    stencil[3] = (i  , j-1, k  )
    stencil[4] = (i  , j  , k-1)
    return stencil, nb
end

# Stencil for bottom south-west boundary cell
function get_stencil!(cell, stencil_type::BottomSouthWestBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 3
    #stencil = Array{Tuple,1}(undef,3)
    stencil[1] = (i+1, j  , k  )
    stencil[2] = (i  , j+1, k  )
    stencil[3] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for bottom south-east boundary cell
function get_stencil!(cell, stencil_type::BottomSouthEastBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 3
    #stencil = Array{Tuple,1}(undef,3)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i  , j+1, k  )
    stencil[3] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for bottom north-west boundary cell
function get_stencil!(cell, stencil_type::BottomNorthWestBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 3
    #stencil = Array{Tuple,1}(undef,3)
    stencil[1] = (i+1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for bottom north-east boundary cell
function get_stencil!(cell, stencil_type::BottomNorthEastBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 3
    #stencil = Array{Tuple,1}(undef,3)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j  , k+1)
    return stencil, nb
end

# Stencil for top south-west boundary cell
function get_stencil!(cell, stencil_type::TopSouthWestBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 3
    #stencil = Array{Tuple,1}(undef,3)
    stencil[1] = (i+1, j  , k  )
    stencil[2] = (i  , j+1, k  )
    stencil[3] = (i  , j  , k-1)
    return stencil, nb
end

# Stencil for top south-east boundary cell
function get_stencil!(cell, stencil_type::TopSouthEastBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 3
    #stencil = Array{Tuple,1}(undef,3)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i  , j+1, k  )
    stencil[3] = (i  , j  , k-1)
    return stencil, nb
end

# Stencil for top north-west boundary cell
function get_stencil!(cell, stencil_type::TopNorthWestBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 3
    #stencil = Array{Tuple,1}(undef,3)
    stencil[1] = (i+1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j  , k-1)
    return stencil, nb
end

# Stencil for top north-east boundary cell
function get_stencil!(cell, stencil_type::TopNorthEastBoundaryCell, stencil, nb)
    (i, j, k) = cell
    nb = 3
    #stencil = Array{Tuple,1}(undef,3)
    stencil[1] = (i-1, j  , k  )
    stencil[2] = (i  , j-1, k  )
    stencil[3] = (i  , j  , k-1)
    return stencil, nb
end




function set_stencil_type(grid) #, stencil_type)

    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    stencil_type = Array{Any,3}(undef, (nx,ny,nz))

    for k = 2:nz
        for j = 2:ny
            for i = 1:nx
                stencil_type[i,j,k] = InternalCell()
            end
        end
    end

    # Initialize all cells to be of Internal cell type
    #stencil_type[: ,: ,: ] .= InternalCell

    # Modify the 6 faces of the domain
    for k = 2:nz
        for j = 2:ny
            stencil_type[1 ,j ,k ] = WestBoundaryCell()
            stencil_type[nx,j ,k ] = EastBoundaryCell()
        end
    end

    for k = 2:nz
        for i = 1:nx
            stencil_type[i ,1 ,k ] = SouthBoundaryCell()
            stencil_type[i ,ny,k ] = NorthBoundaryCell()
        end
    end

    for j = 2:ny
        for i = 1:nx
            stencil_type[i ,j ,1 ] = BottomBoundaryCell()
            stencil_type[i ,j ,nz] = TopBoundaryCell()
        end
    end

    # Modify the 12 sides of the domain
    for k = 1:nz
        stencil_type[1 ,1 ,k ] = SouthWestBoundaryCell()
        stencil_type[nx,1 ,k ] = SouthEastBoundaryCell()
        stencil_type[1 ,ny,k ] = NorthWestBoundaryCell()
        stencil_type[nx,ny,k ] = NorthEastBoundaryCell()
    end

    for j = 1:ny
        stencil_type[1 ,j ,1 ] = BottomWestBoundaryCell()
        stencil_type[nx,j ,1 ] = BottomEastBoundaryCell()
        stencil_type[1 ,j ,nz] = TopWestBoundaryCell()
        stencil_type[nx,j ,nz] = TopEastBoundaryCell()
    end

    for i = 1:nx
        stencil_type[i ,1 ,1 ] = BottomSouthBoundaryCell()
        stencil_type[i ,ny,1 ] = BottomNorthBoundaryCell()
        stencil_type[i ,1 ,nz] = TopSouthBoundaryCell()
        stencil_type[i ,ny,nz] = TopNorthBoundaryCell()
    end

    # Finally modify the 4 domain corners
    stencil_type[1 ,1 ,1 ] = BottomSouthWestBoundaryCell()
    stencil_type[nx,1 ,1 ] = BottomSouthEastBoundaryCell()
    stencil_type[1 ,ny,1 ] = BottomNorthWestBoundaryCell()
    stencil_type[nx,ny,1 ] = BottomNorthEastBoundaryCell()
    stencil_type[1 ,1 ,nz] = TopSouthWestBoundaryCell()
    stencil_type[nx,1 ,nz] = TopSouthEastBoundaryCell()
    stencil_type[1 ,ny,nz] = TopNorthWestBoundaryCell()
    stencil_type[nx,ny,nz] = TopNorthEastBoundaryCell()

    return stencil_type
end

@inline function kernel!(cell, h, hnew, K, stencil, nb)
    (i, j, k ) = cell
    #@show nb
    hnew[i,j,k] = h[i,j,k]
    for ib = 1:nb
        hnew[i,j,k] += add_flux(cell, h, K, stencil[ib])
    end
    return hnew
end

function kernel_no_stencil!(cell, h, hnew, K)
    (i, j, k ) = cell

    hnew[i,j,k] = h[i,j,k] +
                    K[i,j,k]*(h[i+1,j  ,k  ]-2h[i,j,k] + h[i-1,j  ,k  ] +
                              h[i  ,j+1,k  ]-2h[i,j,k] + h[i  ,j-1,k  ] +
                              h[i  ,j  ,k+1]-2h[i,j,k] + h[i  ,j  ,k-1] )
    return hnew
end

function add_flux(cell, h, K, neighb)
    (i ,j ,k ) = cell
    (in,jn,kn) = neighb
    # Take the direction of the neighbour into account
    s = Float64(sign(i-in)+sign(j-jn)+sign(k-kn))
    flux = s * K[i,j,k] * (h[i,j,k]-h[in,jn,kn])
end

function test_full_stencil_based!(grid, h, hnew, K, stencil_type)
    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    nb      = 6
    stencil = Array{Tuple,1}(undef,6)

    @inbounds for k = 1:nz
        @inbounds for j = 1:ny
            @inbounds for i = 1:nx
                cell = (i,j,k)
                stencil, nb = get_stencil!(cell, stencil_type[i,j,k], stencil, nb)
                hnew = kernel!(cell, h, hnew, K, stencil, nb)
                #@show size(neighbours)
            end
        end
    end

end

function test_stencil_based!(grid, h, hnew, K, stencil_type)
    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    nb      = 6
    stencil = Array{Tuple,1}(undef,6)

    @inbounds for k = 2:nz-1
        @inbounds for j = 2:ny-1
            @inbounds for i = 2:nx-1
                cell = (i,j,k)
                stencil, nb = get_stencil!(cell, stencil_type[i,j,k], stencil, nb)
                hnew = kernel!(cell, h, hnew, K, stencil, nb)
            end
        end
    end

end

function test_no_stencil!(grid, h, hnew, K)
    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    @inbounds for k = 2:nz-1
        @inbounds for j = 2:ny-1
            @inbounds for i = 2:nx-1
                cell = (i,j,k)
                hnew = kernel_no_stencil!(cell, h, hnew, K)
            end
        end
    end

end

function test_stencil()

    nx = 100
    ny = 80
    nz = 50

    @show total = nx*ny*nz
    @show internal = (nx-1)*(ny-1)*(nz-1)
    @show boundary = total - internal
    @show boundary_β = 100*(boundary/total)

    dx = 10.0
    dy = 10.0
    dz = 10.0
    useCUDA = false

    h0 = 10.0
    K0 = 1.0e-4

    h = fill(h0,(nx,ny,nz))
    hnew = copy(h)
    K = fill(K0,(nx,ny,nz))

    x,y,z = grid_coords(nx, ny, nz, dx, dy, dz, useCUDA)

    grid = Grid(nx,ny, nz, dx, dy, dz, x, y, z)

    #stencil_type = fill(nothing, (nx, ny, nz))

    stencil_type = set_stencil_type(grid) #, stencil_type)

    @time test_full_stencil_based!(grid, h, hnew, K, stencil_type)
    @time test_stencil_based!(grid, h, hnew, K, stencil_type)
    @time test_no_stencil!(grid, h, hnew, K)

end

test_stencil()
