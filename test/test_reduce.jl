using CUDAnative
using CuArrays

@testset "reduce_max" begin
    N = 160
    u = rand(N, N)
    u_d = CuArray(u)

    ths = (16, 16)
    bls = (Int(ceil(N / ths[1])), Int(ceil(N / ths[2])))
    @cuda threads = ths blocks = bls reduceMax!(u_d)
    @test u_d[1,1] == maximum(u)
end
