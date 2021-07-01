using SparseArrays

mutable struct MySparseMatrix{FT, IntAT, FloatAT}
    A::SparseMatrixCSC{FT, Int64}
    rowptr::IntAT
    colptr::IntAT
    vals::FloatAT
end

"""Initialize the matrix and right hand side used in the linear system for a domain of size (nx,ny,nz)."""
function initLinearSystem(nx, ny, nz, myFloat, useCUDA)
    # The total number of nonzero entries for a 7-diagonal diffusion matrix
    N =     7 * (nx - 2) * (ny - 2) * (nz - 2) + 
        2 * 6 * (nx - 2) * (ny - 2) + 
        2 * 6 * (nx - 2) * (nz - 2) +
        2 * 6 * (ny - 2) * (nz - 2) + 
        2 * 6 * (nx - 2) * (nz - 2) +
        4 * 5 * (nx - 2) + 
        4 * 5 * (ny - 2) + 
        4 * 5 * (nz - 2) +  
        8 * 4

    # Fill the row and column pointers and the value vector with dummy values
    rowptr = [i for i = 1:N]
    colptr = [i for i = 1:N]
    vals    = [myFloat(1.0) for i = 1:N]

    # Generate the sparse matrix from the vectors
    M  = sparse(rowptr, colptr, vals)
    IntAT = typeof(rowptr)
    FloatAT = typeof(vals)
    A = MySparseMatrix{myFloat, IntAT, FloatAT}(M, rowptr, colptr, vals)
    
    # Also make the right-hand-side vector with the size of the total number of equations
    b = Vector{myFloat}(1:nx*ny*nz)

    return A, b
end

"""Initialize the preconditioner."""
function setPreconditioner(thePreconditioner, A)
    if thePreconditioner == "DiagonalPreconditioner"
        # Diagonal preconditioner
        p = DiagonalPreconditioner(A)
    elseif thePreconditioner == "CholeskyPreconditioner"
        # Incomplete Cholesky preconditioner with cut-off level 2
        p = CholeskyPreconditioner(A, 2)
    elseif thePreconditioner == "AMG_RugeStubenPreconditioner"
        # Algebraic multigrid preconditioner (AMG): Ruge-Stuben variant
        p = AMGPreconditioner{RugeStuben}(A)
    elseif thePreconditioner == "AMG_SmoothedAggregationPreconditioner"
        # Algebraic multigrid preconditioner (AMG) using Smoothed aggregation 
        p = AMGPreconditioner{SmoothedAggregation}(A)
    elseif thePreconditioner == "AMG_RugeStuben2Preconditioner"
        # Algebraic multigrid preconditioner (AMG): Ruge-Stuben variant 2
        # (from Preconditioners package)
        p = aspreconditioner(ruge_stuben(A))
    end
    return p
end

"""Update the already existing preconditioner."""
function setPreconditioner!(p, thePreconditioner, A)
    if thePreconditioner == "DiagonalPreconditioner"
        # Diagonal preconditioner
        p = UpdatePreconditioner!(p, A)
    elseif thePreconditioner == "CholeskyPreconditioner"
        # Incomplete Cholesky preconditioner with cut-off level 2
        p = UpdatePreconditioner!(p, A, 2)
    elseif thePreconditioner == "AMG_RugeStubenPreconditioner"
        # Algebraic multigrid preconditioner (AMG): Ruge-Stuben variant
        p = UpdatePreconditioner!(p, A)
    elseif thePreconditioner == "AMG_SmoothedAggregationPreconditioner"
        # Algebraic multigrid preconditioner (AMG) using Smoothed aggregation 
        p = UpdatePreconditioner!(p, A)
    elseif thePreconditioner == "AMG_RugeStuben2Preconditioner"
        # Algebraic multigrid preconditioner (AMG): Ruge-Stuben variant 2
        # (from Preconditioners package)
        p = UpdatePreconditioner!(p,A)
    end
    return p
end