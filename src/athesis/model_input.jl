# model_input.jl
using Adapt

abstract type ModelInput end

struct Model_input_2dh <: ModelInput
    model_type::String
    data_type #::AbstractArray
    ndim::Int64
    nx::Int64
    ny::Int64
    Δx::Float64
    Δy::Float64
    Δt::Float64
    tend::Float64
    K0::Float64
    h0::Float64
    u0::Float64
    v0::Float64
    source::Float64
    i_src::Int64
    j_src::Int64
    duration::Float64
end

struct Model_input_2dv <: ModelInput
    model_type::String
    data_type
    ndim::Int64
    nx::Int64
    nz::Int64
    Δx::Float64
    Δz::Float64
    Δt::Float64
    tend::Float64
    K0::Float64
    h0::Float64
    u0::Float64
    w0::Float64
    source::Float64
    i_src::Int64
    k_src::Int64
    duration::Float64
end

struct Model_input_3d <: ModelInput
    model_type::String
    data_type
    ndim::Int64
    nx::Int64
    ny::Int64
    nz::Int64
    Δx::Float64
    Δy::Float64
    Δz::Float64
    Δt::Float64
    tend::Float64
    K0::Float64
    h0::Float64
    u0::Float64
    v0::Float64
    w0::Float64
    source::Float64
    i_src::Int64
    j_src::Int64
    k_src::Int64
    duration::Float64
end

#
# This function takes the model input
# Creates the model
#
function model_input()

    println("Reading model input ...")

    model_type = "2DH" # in large caps, please!

    # Model size per dimension
    nx   = 200
    ny   = 100
    nz   = 20

    # Grid sizes
    Δx   = 10.0
    Δy   = 10.0
    Δz   = 10.0

    # Time step
    Δt   = 100.0

    # hydraulic_conductivity
    K0   = 1.0e-3

    # Initial condition
    h0   = 10.0
    u0   = 0.0
    v0   = 0.0
    w0   = 0.0

    # Source/discharge data
    i_src  = 100
    j_src  = 50
    k_src  = 5
    duration = 30000.0
    source = 1.0e-2

    # Simulation end time
    tend = 100000.0

    # Backend selection
    println("Hit c for cuda...")

    # For debugging purposes: set c here instead of using readline()
    c = readline()
    #c = "a"
    useCUDA = (c == "c")

    # Set corresponding data type
    if useCUDA
        data_type = CuArray
    else
        data_type = Array
    end

    # Store the input in tuple "input"
    if model_type == "2DH"
        ndim = 2
        println("Grid specified: ", model_type, " grid with\n nx = ", nx, ",\n ny = ", ny, " and\n Δx = ", Δx, ",\n Δy = ", Δy)
        input = Model_input_2dh(model_type, data_type, ndim, nx, ny, Δx, Δy, Δt, tend, K0, h0, u0, v0, source, i_src, j_src, duration)
    elseif model_type == "2DV"
        ndim = 2
        println("Grid specified: ", model_type, " grid with\n nx = ", nx, ",\n nz = ", nz, " and\n Δx = ", Δx, ",\n Δz = ", Δz)
        input = Model_input_2dv(model_type, data_type, ndim, nx, nz, Δx, Δz, Δt, tend, K0, h0, u0, w0, source, i_src, k_src, duration)
    elseif model_type == "3D"
        ndim = 3
        println("Grid specified: ", model_type, " grid with\n nx = ", nx, ",\n ny = ", ny, ",\n nz = ", nz, " and\n Δx = ", Δx, ",\n Δy = ", Δy, ",\n Δz = ", Δz)
        input = Model_input_3d(model_type, data_type, ndim, nx, ny, nz, Δx, Δy, Δz, Δt, tend, K0, h0, u0, v0, w0, source, i_src, j_src, k_src, duration)
    end

end
