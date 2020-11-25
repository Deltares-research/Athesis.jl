abstract type Modeltype
end

# mutable struct myFloat <:AbstractFloat
# end

mutable struct Parameters{AT, FT}
    K::AT
    specificStorage::FT
    #Δt::Float64
end

mutable struct BoundaryConditions{AT}
    bcPressure::AT
end

mutable struct Physics{AT, FT}
    g::FT
    nuₕ::AT
    nuᵥ::AT
end

mutable struct ModelComponents
    #advection::Bool
    #diffusion::Bool
    #bottom_friction::Bool
    #wind_friction::Bool
    #gravity::Bool
    hydraulicConductivity::Bool
    massFlux::Bool
    pressure::Bool
    sources::Bool
end

mutable struct State{AT}
    h::AT
    u::AT
    v::AT
    w::AT
    hⁿ⁺¹::AT
    uⁿ⁺¹::AT
    vⁿ⁺¹::AT
    wⁿ⁺¹::AT
end

mutable struct Source{AT, FT}
    i_src::Int64
    j_src::Int64
    k_src::Int64
    duration::FT
    discharge::FT
    externalSource::AT
end

mutable struct Recharge{FT}
    #duration::AbstractFloat
    constRecharge::FT
    rechargeFlux::FT
    rechargeFactor::FT
end

mutable struct Model
    source::Source
    recharge::Recharge
    boundaryConditions::BoundaryConditions
end

mutable struct TimeData{FT}
    Δt::FT
    tend::FT
    time::FT
    maxsteps::Int64
end

mutable struct SolverData{AT, FT}
    ΔhConv::FT
    Δh::AT
end

mutable struct Simulation
    grid
    model
    state
    parameters
    timeData
    solverData
end
