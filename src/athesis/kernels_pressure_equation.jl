using KernelAbstractions

@inline mat2vec(i, j, k, nx, ny) = i + (j - 1) * nx + (k - 1) * nx*ny

function averageK(K, i, j, k)
    Kw = harmonicMean(K[i - 1, j,     k    ],     K[i,     j,     k    ])
    Ke = harmonicMean(K[i,     j,     k    ],     K[i + 1, j,     k    ])
    Ks = harmonicMean(K[i,     j - 1, k    ],     K[i,     j,     k    ])
    Kn = harmonicMean(K[i,     j,     k    ],     K[i,     j + 1, k    ])
    Kb = harmonicMean(K[i,     j,     k - 1],     K[i,     j,     k    ])
    Kt = harmonicMean(K[i,     j,     k    ],     K[i,     j,     k + 1])
    return Kw, Ke, Ks, Kn, Kb, Kt
end

function harmonicMean(v1, v2)
    return 2.0 * v1 * v2 / (v1 + v2)
end

@kernel function explicitPressureKernel!(source, h, hⁿ⁺¹, Δx, Δy, Δz, Δt, K, SS)
    i, j, k = @index(Global, NTuple)
    
    @inbounds begin
        F = K[i,j,k] * (
            (h[i + 1,j,k] + h[i - 1,j,k] - 2.0 * h[i,j,k]) / (Δx * Δx) +
            (h[i,j + 1,k] + h[i,j - 1,k] - 2.0 * h[i,j,k]) / (Δy * Δy) +
            (h[i,j,k + 1] + h[i,j,k - 1] - 2.0 * h[i,j,k]) / (Δz * Δz)
            ) + source[i,j,k] / (Δx * Δy)

        hⁿ⁺¹[i,j,k] = h[i,j,k] + Δt * F * (1.0 / SS)
    end

end

@kernel function implicitPressureKernel!(myA::MySparseMatrix, rhs, source, h, hⁿ⁺¹, nx, ny, nz, Δx, Δy, Δz, Δt, K, SS)
    i, j, k = @index(Global, NTuple)

    Δx² = Δx * Δx
    Δy² = Δy * Δy
    Δz² = Δz * Δz

    @inbounds begin

        # Present cell/equation index
        ijk = mat2vec(i, j, k, nx, ny)

        # Lower diagonal from negative z-dir
        if k > 1
            ijk1 = mat2vec(i, j , k - 1, nx, ny)
            myA.rowptr[ijk] = ijk
            myA.colptr[ijk] = ijk1
            myA.vals[ijk]   = -K[i,j,k] / Δz²
        end

        # Lower diagonal from negative y-dir
        if j > 1
            ijk1 = mat2vec(i, j - 1, k, nx, ny)
            myA.rowptr[ijk] = ijk
            myA.colptr[ijk] = ijk1
            myA.vals[ijk]   = -K[i,j,k] / Δy²
        end

        # Lower diagonal from negative x-dir
        if i > 1
            ijk1 = mat2vec(i - 1, j, k, nx, ny)
            myA.rowptr[ijk] = ijk
            myA.colptr[ijk] = ijk1
            myA.vals[ijk]   = -K[i,j,k] / Δx²
        end

        # The main diagonal
        myA.rowptr[ijk] = ijk
        myA.colptr[ijk] = ijk
        myA.vals[ijk]   = 2.0 * K[i,j,k] * (1.0 / Δx² + 1.0 / Δy² + 1.0 / Δz²)

        # Upper diagonal from positive x-dir
        if i < nx
            ijk1 = mat2vec(i + 1, j, k, nx, ny)
            myA.rowptr[ijk] = ijk
            myA.colptr[ijk] = ijk1
            myA.vals[ijk]    = -K[i,j,k] / Δx²
        end

        # Upper diagonal from positive y-dir
        if j < ny
            ijk1 = mat2vec(i, j + 1, k, nx, ny)
            myA.rowptr[ijk] = ijk
            myA.colptr[ijk] = ijk1
            myA.vals[ijk]   = -K[i,j,k] / Δy²
        end

        # Upper diagonal from positive z-dir
        if j < nz
            ijk1 = mat2vec(i, j, k + 1, nx, ny)
            myA.rowptr[ijk] = ijk
            myA.colptr[ijk] = ijk1
            myA.vals[ijk]   = -K[i,j,k] / Δz²
        end
    end

end

@inline function pressureKernelAveraged!(i, j, k, source, h, u, v, w, hⁿ⁺¹, uⁿ⁺¹, vⁿ⁺¹, wⁿ⁺¹, Δx, Δy, Δz, Δt, K, SS)

    Kw, Ke, Ks, Kn, Kb, Kt = averageK(K, i, j, k)

    F = (1.0 / Δx * Δx) * (Ke * (h[i + 1,j  ,k  ] - h[i,j,k]) - Kw * (h[i,j,k] - h[i - 1,j  ,k  ])) +
        (1.0 / Δy * Δy) * (Kn * (h[i  ,j + 1,k  ] - h[i,j,k]) - Ks * (h[i,j,k] - h[i  ,j - 1,k  ])) +
        (1.0 / Δz * Δz) * (Kt * (h[i  ,j  ,k + 1] - h[i,j,k]) - Kb * (h[i,j,k] - h[i  ,j  ,k - 1])) +
        source[i,j,k] / (Δx * Δy)

    hⁿ⁺¹[i,j,k] = h[i,j,k] + Δt * F * (1.0 / SS)
end
