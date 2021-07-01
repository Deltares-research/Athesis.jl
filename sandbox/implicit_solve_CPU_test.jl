
# solve a linear system of equations
using Plots
using OffsetArrays
using CUDA

# using Krylov
using LinearAlgebra
using IterativeSolvers
using SparseArrays
using Preconditioners
using TimerOutputs
# using AlgebraicMultigrid

macro synctimeit(timer, label, body)
    if has_cuda()
        return :( 
            if ($(esc(timer)).enabled) 
                @timeit $(esc(timer)) $(esc(label)) CUDA.@sync $(esc(body)) 
            else
                $(esc(body))
            end
            )
    else
        return :( 
            if ($(esc(timer)).enabled) 
                @timeit $(esc(timer)) $(esc(label)) $(esc(body)) 
            else
                $(esc(body))
            end
            )
    end

end

# abstract type ThePreconditioner end
# mutable struct DiagPreconditioner <: ThePreconditioner end
# mutable struct CholeskyPreconditioner <: ThePreconditioner end
# mutable struct AMG_RugeStubenPreconditioner <: ThePreconditioner end
# mutable struct AMG_SmoothedAggregationPreconditioner <: ThePreconditioner end
# mutable struct AMG_RugeStuben2Preconditioner <: ThePreconditioner end

mutable struct MySparseMatrix{FT, IntAT, FloatAT}
    A::SparseMatrixCSC{FT, Int64}
    rowptr::IntAT
    colptr::IntAT
    vals::FloatAT
end

@inline mat2vec(i, j, nx) = i + (j - 1) * nx

"""Create a sparse diffusion matrix in Compressed Column Storage for the CPU"""
function diffusion_matrix_CSC_CPU(nx, ny, dx, dy) # A::SparseMatrixCSC{Int64, Int64},
    # For each equation/row in the matrix we fill 5 entries:
    # The main diagonal and 4 off-diagonals
    # Start with all internal points

    dx2 = dx * dx
    dy2 = dy * dy

    N = 5 * (nx - 2) * (ny - 2) + 4 * 2 * (nx - 2) + 4 * 2 * (ny - 2) + 3 * 4
    #rowptr = zeros(Int64, N)
    #colptr = zeros(Int64, N)
    #val    = zeros(Float64, N)

    # rowptr = CuArray(rowptr)
    # colptr = CuArray(colptr)
    # val    = CuArray(val)

    idx = 0
    for j = 1:ny
        for i = 1:nx
            ij = mat2vec(i, j, nx)

            # Lower diagonal from negative y-dir
            if j > 1
                # ij1 = mat2vec(i,j-1,nx)
                # A[ij,ij1] = -1.0/dy2
                idx += 1
                rowvals(A)[idx] = ij
                nzrange(A)[idx] = mat2vec(i, j - 1, nx)
                nonzeros(A)[idx] = -1.0 / dy2
            end

            # Lower diagonal from negative x-dir
            if i > 1
                # ij1 = mat2vec(i-1,j,nx)
                # A[ij,ij1] = -1.0/dx2
                idx += 1
                rowvals(A)[idx] = ij
                nzrange(A)[idx] = mat2vec(i - 1, j, nx)
                nonzeros(A)[idx] = -1.0 / dx2
            end

            # The main diagonal
            # A[ij,ij] = 2.0/dx2 + 2.0/dy2
            idx += 1
            rowvals(A)[idx] = ij
            nzrange(A)[idx] = ij
            nonzeros(A)[idx] = 2.0 / dx2 + 2.0 / dy2

            # Upper diagonal from positive x-dir
            if i < nx
                # ij1 = mat2vec(i+1,j,nx)
                # A[ij,ij1] = -1.0/dx2
                idx += 1
                rowvals(A)[idx] = ij
                nzrange(A)[idx] = mat2vec(i + 1, j, nx)
                nonzeros(A)[idx] = -1.0 / dx2
            end

            # Upper diagonal from positive y-dir
            if j < ny
                # ij1 = mat2vec(i,j+1,nx)
                # A[ij,ij1] = -1.0/dy2
                idx += 1
                rowvals(A)[idx] = ij
                nzrange(A)[idx] = mat2vec(i, j + 1, nx)
                nonzeros(A)[idx] = -1.0 / dy2
            end
            # @show idx
        end
    end

    # return sparse CSC matrix
    return A # = sparse(rowptr, colptr, val)

end

# @kernel function kernel_diffusion_matrix!(myA::)
    
function diffusion_matrix_CSC!(myA::MySparseMatrix, nx, ny, dx, dy)
    # For each equation/row in the matrix we fill 5 entries:
    # The main diagonal and 4 off-diagonals
    # Start with all internal points

    dx2 = dx * dx
    dy2 = dy * dy

    #N = 5 * (nx - 2) * (ny - 2) + 4 * 2 * (nx - 2) + 4 * 2 * (ny - 2) + 3 * 4
    #rowptr = zeros(Int64, N)
    #colptr = zeros(Int64, N)
    #val    = zeros(Float64, N)

    idx = 0
    for j = 1:ny
        for i = 1:nx
            ij = mat2vec(i, j, nx)

            # Lower diagonal from negative y-dir
            if j > 1
                ij1 = mat2vec(i, j - 1, nx)
                # A[ij,ij1] = -1.0/dy2
                idx += 1
                myA.rowptr[idx] = ij
                myA.colptr[idx] = ij1
                myA.vals[idx]   = -1.0 / dy2
            end

            # Lower diagonal from negative x-dir
            if i > 1
                ij1 = mat2vec(i - 1, j, nx)
                # A[ij,ij1] = -1.0/dx2
                idx += 1
                myA.rowptr[idx] = ij
                myA.colptr[idx] = ij1
                myA.vals[idx]   = -1.0 / dx2
            end

            # The main diagonal
            # A[ij,ij] = 2.0/dx2 + 2.0/dy2
            idx += 1
            myA.rowptr[idx] = ij
            myA.colptr[idx] = ij
            myA.vals[idx]   = 2.0 / dx2 + 2.0 / dy2

            # Upper diagonal from positive x-dir
            if i < nx
                ij1 = mat2vec(i + 1, j, nx)
                # A[ij,ij1] = -1.0/dx2
                idx += 1
                myA.rowptr[idx] = ij
                myA.colptr[idx] = ij1
                myA.vals[idx]    = -1.0 / dx2
            end

            # Upper diagonal from positive y-dir
            if j < ny
                ij1 = mat2vec(i, j + 1, nx)
                # A[ij,ij1] = -1.0/dy2
                idx += 1
                myA.rowptr[idx] = ij
                myA.colptr[idx] = ij1
                myA.vals[idx]   = -1.0 / dy2
            end
            # @show idx
        end
    end

    # return sparse CSC matrix
    return myA.A = sparse(myA.rowptr, myA.colptr, myA.vals)

end

function set_rhs!(n, b, nx, ny)
    for j = 1:ny
        for i = 1:nx
            # b[i] = Float64(i)/Float64(length(b))
            if i == Int64(nx / 2) && j == Int64(ny / 2)
                ij = mat2vec(i, j, nx)
                b[ij] = n
            end
        end
    end
    return b
end

function plot_solution(x, nx, ny)
    y = reshape(x, (nx, ny))
    p1 = contour(y, fill=true)
    display(p1)
end

function init_sparse_matrix(nx, ny)
    N = 5 * (nx - 2) * (ny - 2) + 4 * 2 * (nx - 2) + 4 * 2 * (ny - 2) + 3 * 4
    rowptr = [i for i = 1:N]
    colptr = [i for i = 1:N]
    vals    = zeros(Float64, N)

    #rowptr = [1,1]
    #colptr = [1,1]
    #val    = [0.0, 0.0]
    M  = sparse(rowptr, colptr, vals)
    FT = typeof(vals[1]) 
    IntAT = typeof(rowptr)
    FloatAT = typeof(vals)
    A = MySparseMatrix{FT, IntAT, FloatAT}(M, rowptr, colptr, vals)
end

function setPreconditioner(thePreconditioner, A)
    if typeof(thePreconditioner) == DiagonalPreconditioner
        # Diagonal preconditioner
        p = DiagonalPreconditioner(A)
    elseif typeof(thePreconditioner) == CholeskyPreconditioner
        # Incomplete Cholesky preconditioner with cut-off level 2
        p = CholeskyPreconditioner(A, 2)
    elseif typeof(thePreconditioner) == AMG_RugeStubenPreconditioner
        # Algebraic multigrid preconditioner (AMG): Ruge-Stuben variant
        p = AMGPreconditioner{RugeStuben}(A)
    elseif typeof(thePreconditioner) == AMG_SmoothedAggregationPreconditioner
        # Algebraic multigrid preconditioner (AMG) using Smoothed aggregation 
        p = AMGPreconditioner{SmoothedAggregation}(A)
    elseif typeof(thePreconditioner) == AMG_RugeStuben2Preconditioner
        # Algebraic multigrid preconditioner (AMG): Ruge-Stuben variant 2
        # (from Preconditioners package)
        p = aspreconditioner(ruge_stuben(A))
    end
    return p
end

function setPreconditioner!(p, thePreconditioner, A)
    if typeof(thePreconditioner) == DiagonalPreconditioner
        # Diagonal preconditioner
        return p = UpdatePreconditioner!(p, A)
    elseif typeof(thePreconditioner) == CholeskyPreconditioner
        # Incomplete Cholesky preconditioner with cut-off level 2
        return p = UpdatePreconditioner!(p, A, 2)
    elseif typeof(thePreconditioner) == AMG_RugeStubenPreconditioner
        # Algebraic multigrid preconditioner (AMG): Ruge-Stuben variant
        return p = UpdatePreconditioner!(p, A)
    elseif typeof(thePreconditioner) == AMG_SmoothedAggregationPreconditioner
        # Algebraic multigrid preconditioner (AMG) using Smoothed aggregation 
        return p = UpdatePreconditioner!(p, A)
    elseif typeof(thePreconditioner) == AMG_RugeStuben2Preconditioner
        # Algebraic multigrid preconditioner (AMG): Ruge-Stuben variant 2
        # (from Preconditioners package)
        return p = UpdatePreconditioner!(p,A)
    end
end

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

#"""Function to test iterative solvers and preconditioners on CPU"""
function test_solver()
    
    # Domain size
    nx = 500
    ny = 500
    dx = 1.0
    dy = 1.0

    numsteps = 10
    tol = 1.0e-8

    first = true

    # Matrix size
    @show N = nx * ny

    # Solution vector
    x = zeros(N)

    # Right hand side
    b = zeros(N)
    
    # Initialize sparse matrix
    myA = init_sparse_matrix(nx,ny)
    A = diffusion_matrix_CSC!(myA, nx, ny, dx, dy)

    # Set and initialize the preconditioner
    #myPreconditioner = AMG_SmoothedAggregationPreconditioner()
    myPreconditioner = "AMG_SmoothedAggregationPreconditioner"
    p = setPreconditioner(myPreconditioner, A)

    to = TimerOutput()

    # Now perform a number of (time) steps
    for n = 1:numsteps
        # Set a specific right hand side
        b = set_rhs!(n, b, nx, ny)

        # Now form a sparse matrix from diagonals
        @synctimeit to "matrix setup" begin
            A = diffusion_matrix_CSC!(myA, nx, ny, dx, dy)
        end

        # Set up the preconditioner
        @synctimeit to "preconditioner" begin
            p = setPreconditioner!(p, myPreconditioner, A)
        end

        # Solve the system of equations
        @synctimeit to "solve" begin
            x, ch = cg!(x, A, b, Pl=p, abstol=tol, maxiter=1000, verbose=false, log=true)
            @show ch
            @show ch[:resnorm, ch.iters]

            #fig = plot!(ch, :resnorm, sep = :blue)
            #display(fig)
        end

        # Plot result
        #plot_solution(x, nx, ny)

    end

    # Show timer output
    print_timer(to)
    
    
end

@time test_solver()