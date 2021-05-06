using KernelAbstractions

function averageK(K, i, j, k)
    Kw = harmonicMean(K[i-1,j  ,k  ],K[i  ,j  ,k  ])
    Ke = harmonicMean(K[i  ,j  ,k  ],K[i+1,j  ,k  ])
    Ks = harmonicMean(K[i  ,j-1,k  ],K[i  ,j  ,k  ])
    Kn = harmonicMean(K[i  ,j  ,k  ],K[i  ,j+1,k  ])
    Kb = harmonicMean(K[i  ,j  ,k-1],K[i  ,j  ,k  ])
    Kt = harmonicMean(K[i  ,j  ,k  ],K[i  ,j  ,k+1])
    return Kw, Ke, Ks, Kn, Kb, Kt
end

function harmonicMean(v1, v2)
    return 2.0*v1*v2/(v1+v2)
end

@kernel function pressureKernel!(source, h, hⁿ⁺¹, Δx, Δy, Δz, Δt, K, SS)
    i, j, k = @index(Global, NTuple)
    
    @inbounds begin
        F = K[i,j,k] * (
            (h[i+1,j,k] + h[i-1,j,k] - 2.0*h[i,j,k])/(Δx*Δx) +
            (h[i,j+1,k] + h[i,j-1,k] - 2.0*h[i,j,k])/(Δy*Δy) +
            (h[i,j,k+1] + h[i,j,k-1] - 2.0*h[i,j,k])/(Δz*Δz)
            ) + source[i,j,k]/(Δx*Δy)

        hⁿ⁺¹[i,j,k] = h[i,j,k] + Δt*F*(1.0/SS)
    end

end

@inline function pressureKernelAveraged!(i, j, k, source, h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹, Δx, Δy, Δz, Δt, K, SS)

    Kw, Ke, Ks, Kn, Kb, Kt = averageK(K,i, j, k)

    F = (1.0/Δx*Δx) * (Ke * (h[i+1,j  ,k  ] - h[i,j,k]) - Kw * (h[i,j,k] - h[i-1,j  ,k  ])) +
        (1.0/Δy*Δy) * (Kn * (h[i  ,j+1,k  ] - h[i,j,k]) - Ks * (h[i,j,k] - h[i  ,j-1,k  ])) +
        (1.0/Δz*Δz) * (Kt * (h[i  ,j  ,k+1] - h[i,j,k]) - Kb * (h[i,j,k] - h[i  ,j  ,k-1])) +
        source[i,j,k]/(Δx*Δy)

    hⁿ⁺¹[i,j,k] = h[i,j,k] + Δt*F*(1.0/SS)
end
