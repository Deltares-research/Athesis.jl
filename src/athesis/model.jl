# model.jl

abstract type Modeltype
end


mutable struct Parameters #{T}
    K::AbstractArray{Float64}
    #Δt::Float64
end

mutable struct Physics #{T}
    g::Float64
    nuₕ::AbstractArray{Float64}
    nuᵥ::AbstractArray{Float64}
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

mutable struct State2dh{T}
    h::T
    u::T
    v::T
    hⁿ⁺¹::T
    uⁿ⁺¹::T
    vⁿ⁺¹::T
end

mutable struct State2dv{T}
    h::T
    u::T
    w::T
    hⁿ⁺¹::T
    uⁿ⁺¹::T
    wⁿ⁺¹::T
end

mutable struct State3d{T}
    h::T
    u::T
    v::T
    w::T
    hⁿ⁺¹::T
    uⁿ⁺¹::T
    vⁿ⁺¹::T
    wⁿ⁺¹::T
end

abstract type Source end

mutable struct Source2dh <: Source # {T}
    i_src::Int64
    j_src::Int64
    duration::Float64
    discharge::Float64
    external_source::AbstractArray{Float64,2}
end

mutable struct Source2dv <: Source # {T}
    i_src::Int64
    k_src::Int64
    duration::Float64
    discharge::Float64
    external_source::AbstractArray{Float64,2}
end

mutable struct Source3d <: Source # {T}
    i_src::Int64
    j_src::Int64
    k_src::Int64
    duration::Float64
    discharge::Float64
    external_source::AbstractArray{Float64,3}
end

mutable struct Model #{T}
    #modeltype::Modeltype
    model_type::String
    #state::State
    source::Source
end

mutable struct Time_data
    Δt::Float64
    tend::Float64
    time::Float64
    maxsteps::Int64
    #Courant::Float64
    #time_integration::String
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
