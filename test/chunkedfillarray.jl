using DiskArrays, DiskArrayTools

@testset "Chunked Fill Array" begin
    a = ChunkedFillArray(1, (100,100),(10,10))
    @test eltype(a) == Int
    @test all(==(1), a)
    @test eachchunk(a) isa  DiskArrays.GridChunks
    afloat = ChunkedFillArray{Float64,2}(1,(100,100), (10,10))
    @test eltype(afloat) == Float64
    @test all(==(1.), a)
    @test eachchunk(a) isa DiskArrays.GridChunks
    afloat = ChunkedFillArray{Float64}(1,(100,100), (10,10))
    @test eltype(afloat) == Float64
    @test all(==(1.), a)
    @test eachchunk(a) isa DiskArrays.GridChunks
end

@testset "Chunked Fill Array Chunks provided" begin
    chunks = DiskArrays.GridChunks(DiskArrays.IrregularChunks([0,20,50,100]), DiskArrays.RegularChunks(10,0,100))
    a = ChunkedFillArray(1, (100,100),chunks)
    @test eltype(a) == Int
    @test all(==(1), a)
    @test eachchunk(a) isa  DiskArrays.GridChunks
    afloat = ChunkedFillArray{Float64,2}(1,(100,100), chunks)
    @test eltype(afloat) == Float64
    @test all(==(1.), a)
    @test eachchunk(a) isa DiskArrays.GridChunks
    afloat = ChunkedFillArray{Float64}(1,(100,100), chunks)
    @test eltype(afloat) == Float64
    @test all(==(1.), a)
    @test eachchunk(a) isa DiskArrays.GridChunks
end