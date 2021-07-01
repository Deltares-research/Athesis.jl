# model_input.jl

mutable struct ModelInput{AT,FT}
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

"""Returns default model input"""
function getDefaultInput(myFloat, config_file)
    config = TOML.parsefile(config_file)

    Δx   = config["Lx"] / (config["nx"] + 1)
    Δy   = config["Ly"] / (config["ny"] + 1)
    Δz   = config["Lz"] / (config["nz"] + 1)

    constRecharge = config["unit_recharge_flux"] * Δx * Δy / (24 * 3600) # (m3/s)

    boundaryPressure = [myFloat(config["h_bc_west"]), myFloat(config["h_bc_east"])]

    # Store the input in tuple "input"
    AT = typeof(boundaryPressure)
    input = ModelInput{AT,myFloat}(config["nx"], config["nz"], config["nz"],
                                   config["Lx"], config["Ly"], config["Ly"],
                                   Δx, Δy, Δz,
                                   config["delta_t"], config["tend"],
                                   config["K0"], config["S0"],
                                   config["h0"], config["u0"], config["v0"], config["w0"],
                                   config["source"], config["i_src"], config["j_src"], config["k_src"], config["duration"],
                                   config["delta_h_conv"],
                                   constRecharge, config["recharge_factor"],
                                   boundaryPressure)

end
