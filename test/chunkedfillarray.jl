using DiskArrays, DiskArrayTools

@testset "Chunked Fill Array" begin
    a = ChunkedFillArray(1, (100,100),(10,10))
    @test eltype(a) == Int
    @test all(==(1), a)
    @test eachchunk(a) isa  DiskArrays.GridChunks
end