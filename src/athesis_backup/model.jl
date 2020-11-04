# model.jl

abstract type Modeltype
end


mutable struct Parameters{T}
    K::T
    specific_storage::Float64
    #Δt::Float64
end

mutable struct BoundaryConditions{T}
    bc_pressure::T
end

mutable struct Physics{T}
    g::Float64
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
    duration::Float64
    discharge::Float64
    external_source::T
end

mutable struct Recharge
    #duration::AbstractFloat
    const_recharge::Float64
    recharge_flux::Float64
    recharge_factor::Float64
end

mutable struct Model #{T}
    #modeltype::Modeltype
    #state::State
    source::Source
    recharge::Recharge
    boundary_conditions::BoundaryConditions
end

mutable struct Time_data
    Δt::Float64
    tend::Float64
    time::Float64
    maxsteps::Int64
    #Courant::Float64
    #time_integration::String
end

mutable struct Solver_data{T}
    hclose::Float64
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
    conv_limit::Float64
end
