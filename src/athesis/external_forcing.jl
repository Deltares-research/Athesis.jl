
function setSources!(time, source)
    # 3D implementation
    i_src = source.i_src
    j_src = source.j_src
    k_src = source.k_src
    duration = source.duration
    source.externalSource .= 0.0
    if time < duration
        # Set a single source of at the prescribed point
        source.externalSource[i_src, j_src, k_src] = source.discharge
    else
        source.externalSource[i_src, j_src, k_src] = 0.0
    end
end


function setRecharge!(time, recharge)
    #duration = recharge.duration
    #if time < duration
        # Set a single source of at the prescribed point
        recharge.rechargeFlux = recharge.rechargeFactor*recharge.constRecharge
    #else
    #    recharge.recharge_flux = 0.0
    #end
end
