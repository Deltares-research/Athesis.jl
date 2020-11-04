# parallel reduction: put max at u[1,1]
function reduceMax!(u)
    tx = blockDim().x * (blockIdx().x - 1) + threadIdx().x;
    ty = blockDim().y * (blockIdx().y - 1) + threadIdx().y;

    if tx < size(u,1) + 1 && ty < size(u,2) + 1
        # reduce over x
        stride = blockDim().x >> 1
        while stride > 0
            sync_threads()
            if threadIdx().x <= stride && tx + stride <= size(u,1)
                if u[tx,ty] < u[tx+stride,ty]
                    u[tx,ty] = u[tx+stride,ty]
                end
            end
            stride >>= 1
        end

        # reduce over y
        stride = blockDim().y >> 1
        while stride > 0
            sync_threads()
            if threadIdx().y <= stride && ty + stride <= size(u,2)
                if u[tx,ty] < u[tx,ty+stride]
                    u[tx,ty] = u[tx,ty+stride]
                end
            end
            stride >>= 1
        end
    end

    return nothing
end

# for develop, todo: how to unit test this stuff?
function testReduceMax()
    N = 160
    u = rand(N,N)
    u_d = CuArray(u)

    ths = (16,16)
    bls = (Int(ceil(N/ths[1])),Int(ceil(N/ths[2])))
    @cuda threads=ths blocks=bls reduceMax!(u_d)
    @test u_d[1,1] == maximum(u)
end
