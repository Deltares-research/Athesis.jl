abstract type Modeltype
end


mutable struct Parameters{T}
    K::T
    specificStorage::AbstractFloat
    #Δt::Float64
end

mutable struct BoundaryConditions{T}
    bcPressure::T
end

mutable struct Physics{T}
    g::AbstractFloat
    nuₕ::T
    nuᵥ::T
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

mutable struct State{T}
    h::T
    u::T
    v::T
    w::T
    hⁿ⁺¹::T
    uⁿ⁺¹::T
    vⁿ⁺¹::T
    wⁿ⁺¹::T
end

mutable struct Source{T}
    i_src::Int64
    j_src::Int64
    k_src::Int64
    duration::AbstractFloat
    discharge::AbstractFloat
    externalSource::T
end

mutable struct Recharge
    #duration::AbstractFloat
    constRecharge::AbstractFloat
    rechargeFlux::AbstractFloat
    rechargeFactor::AbstractFloat
end

mutable struct Model
    source::Source
    recharge::Recharge
    boundaryConditions::BoundaryConditions
end

mutable struct TimeData
    Δt::AbstractFloat
    tend::AbstractFloat
    time::AbstractFloat
    maxsteps::Int64
end

mutable struct SolverData{T}
    ΔhConv::AbstractFloat
    Δh::T
end

mutable struct Simulation
    grid
    model
    state
    parameters
    timeData
    solverData
end
