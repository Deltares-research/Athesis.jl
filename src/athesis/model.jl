# model.jl

abstract type Modeltype
end


mutable struct Parameters{T}
    K::T
    specific_storage::AbstractFloat
    #Δt::Float64
end

mutable struct BoundaryConditions{T}
    bc_pressure::T
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
    hydraulic_conductivity::Bool
    massflux::Bool
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
    external_source::T
end

mutable struct Recharge
    #duration::AbstractFloat
    const_recharge::AbstractFloat
    recharge_flux::AbstractFloat
    recharge_factor::AbstractFloat
end

mutable struct Model #{T}
    #modeltype::Modeltype
    #state::State
    source::Source
    recharge::Recharge
    boundary_conditions::BoundaryConditions
end

mutable struct Time_data
    Δt::AbstractFloat
    tend::AbstractFloat
    time::AbstractFloat
    maxsteps::Int64
    #Courant::Float64
    #time_integration::String
end

mutable struct Solver_data{T}
    hclose::AbstractFloat
    Δh::T
end


mutable struct Backend
    GPU::Bool

end
mutable struct TimeIntegration
    name::String
end

mutable struct Solver
    name::String
    conv_limit::AbstractFloat
end
