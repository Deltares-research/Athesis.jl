using OffsetArrays

mutable struct BCTypeDirichlet end

mutable struct BCTypeNeumann end

mutable struct BCTypeRobbins end

mutable struct WaterLevelBC end

mutable struct PressureBC end

mutable struct VelocityBC end

mutable struct FluxBC end

mutable struct BCLocation{T}
    N::Int64
    dim::Int64
    cells::T
end

mutable struct BCValues{T}
N::Int64
    dim::Int64
    values::T
 end

mutable struct BoundaryCondition{string,BCType,BCParameter,BCLocation,BCValues}
    name::string
    type::BCType
    parameter::BCParameter
    location::BCLocation
    values::BCValues
end

ndims = 3
nx  = 10
ny  = 5
nz = 20

dirichlet_bc = BCTypeDirichlet()
neumann_bc   = BCTypeNeumann()
robbins_bc   = BCTypeRobbins()

pressure_bc    = PressureBC()
water_level_bc = WaterLevelBC()
velocity_bc    = VelocityBC()
flux_bc        = FluxBC()

bc_cells  = OffsetArray{Int64}(undef, 0:nx + 1, 0:ny + 1, 0:nz + 1)
bc_values = OffsetArray{Float64}(undef, 0:nx + 1, 0:ny + 1, 0:nz + 1)

bc_cells .= 0
bc_cells .= 0.0

for k = 1:nz
    for j = 1:ny
        bc_cells[0, j,k] = 1
        bc_values[0, j,k] = 1.0
    end
end

location_bc = BCLocation(ny, ndims, bc_cells)
values_bc  = BCValues(ny, ndims, bc_values)

bc_name = "west_bc"

bc_west = BoundaryCondition(bc_name, dirichlet_bc, pressure_bc, location_bc, values_bc)

bc_west.type
bc_west.parameter