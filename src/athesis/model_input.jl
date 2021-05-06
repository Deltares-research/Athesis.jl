# model_input.jl

mutable struct ModelInput{AT, FT}
    nx::Int64
    ny::Int64
    nz::Int64
    Lx::FT
    Ly::FT
    Lz::FT
    Δx::FT
    Δy::FT
    Δz::FT
    Δt::FT
    tend::FT
    K0::FT
    S0::FT
    h0::FT
    u0::FT
    v0::FT
    w0::FT
    source::FT
    i_src::Int64
    j_src::Int64
    k_src::Int64
    duration::FT
    ΔhConv::FT
    constRecharge::FT
    rechargeFactor::FT
    boundaryPressure::AT
end

# returns default model input
function getDefaultInput(myFloat)

    # Model size per dimension (-)
    nx   = 100
    ny   = 100
    nz   = 2

    # Grid extent
    Lx = 80.0
    Ly = 80.0
    Lz = 10.0

    # Grid sizes (m)
    Δx   = Lx/(nx+1)
    Δy   = Ly/(ny+1)
    Δz   = Lz/(nz+1)

    # hydraulic_conductivity (m/s, ~sand)
    K0   = 10.0 / (24*3600)

    # specific storage (1/m)
    S0 = 0.0001

    # Time step (s)
    Δt = S0*(min(Δx,Δy,Δz))^2/(4.0*K0)

    # Initial condition
    h0   = 0.0     # (m)
    u0   = 0.0     # (m/s)
    v0   = 0.0     # (m/s)
    w0   = 0.0     # (m/s)

    # Boundary conditions
    hBCWest = myFloat(1.0)
    hBCEast = myFloat(1.0)
    boundaryPressure = [hBCWest, hBCEast]

    # Source (well) data
    i_src  = 1
    j_src  = 1
    k_src  = 1
    duration = 0.0  # (s)
    source = 0.0    # (m3/s)

    # Recharge data
    # QR_nb = I_nb * M_nb * A_n
    # with QR_nb in (L^3/T, or m3/s)
    # I is the unit rescharge flux (m/s)
    # A is the cell area for cell n
    # M is a multiplier/factor
    constRecharge = 5.0e-4 * Δx * Δy / (24*3600) # (m3/s)
    rechargeFactor = 1.0      # (-)

    # Simulation end time (s)
    tend = 100000.0

    # Convergence criterion for steady state
    ΔhConv = 1e-5

    # Store the input in tuple "input"
    AT = typeof(boundaryPressure)
    input = ModelInput{AT,myFloat}(nx, ny, nz,
                                   Lx, Ly, Lz,
                                   Δx, Δy, Δz,
                                   Δt, tend,
                                   K0, S0,
                                   h0, u0, v0, w0,
                                   source, i_src, j_src, k_src, duration,
                                   ΔhConv,
                                   constRecharge, rechargeFactor,
                                   boundaryPressure)

end
