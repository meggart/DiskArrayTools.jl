using DiskArrays, DiskArrayTools
@testset "ConcatDiskArrays" begin
@testset "1D concatenation" begin
a = [rand(i,3,4) for i = 1:10]
aar = cat(a...,dims=1)
ad = _DiskArray.(a)
aconc = DiskArrayTools.ConcatDiskArray(ad)

@test aconc isa ConcatDiskArray
@test getproperty.(aconc.parents,:parent) == reshape(a,10,1,1)
@test aconc.startinds == ([1, 2, 4, 7, 11, 16, 22, 29, 37, 46],[1],[1])
@test aconc.size == (55,3,4)
@test size(aconc) == (55,3,4)

@test aconc[:,1:2,1:2] == aar[:,1:2,1:2]
@test ad[1].getindex_list == [(1:1,1:2,1:2)]
@test ad[10].getindex_list == [(1:10,1:2,1:2)]
end

@testset "2D concatenation" begin
a = [rand(i,j,4) for i = 1:10, j=1:5]
a1 = [cat(x..., dims = 1) for x in eachcol(a)]
aar = cat(a1...,dims = 2)
ad = _DiskArray.(a)
aconc = ConcatDiskArray(ad)
@test aconc isa ConcatDiskArray
@test getproperty.(aconc.parents,:parent) == reshape(a,10,5,1)
@test aconc.startinds == ([1, 2, 4, 7, 11, 16, 22, 29, 37, 46],[1,2,4,7,11],[1])
@test aconc.size == (55,15,4)
@test size(aconc) == (55,15,4)
@test aconc[:,:,1:2] == aar[:,:,1:2]
@test ad[2,2].getindex_list == [(1:2,1:2,1:2)]
@test ad[8,5].getindex_list == [(1:8,1:5,1:2)]
@test aconc[13:18,3:7,3] == aar[13:18,3:7,3]
end

@testset "Non-matching arrays" begin
    a = Matrix{Matrix{Float64}}(undef,2,2)
    a[1,1] = zeros(2,2)
    a[1,2] = zeros(2,3)
    a[2,1] = zeros(2,3)
    a[2,2] = zeros(2,3)
    @test_throws Exception ConcatDiskArray(_DiskArray.(a))
end

@testset "Operations on chunks" begin
    a = [rand(i,j,4) for i = 1:10, j=1:5]
    a1 = [cat(x..., dims = 1) for x in eachcol(a)]
    aar = cat(a1...,dims = 2)
    ad = _DiskArray.(a)
    aconc = ConcatDiskArray(ad)
    @test isapprox(sum(aconc),sum(aar))
    @test_broken copy(aconc .+ ones(1,15)) == aar .+ ones(1,15)
    @test length(ad[3,5].getindex_list) == 1
end

@testset "Irregular chunks" begin
    a = [rand(i,j,4) for i = 1:10, j=1:5]
    a1 = [cat(x..., dims = 1) for x in eachcol(a)]
    aar = cat(a1...,dims = 2)
    ad = _DiskArray.(a, chunksize = (size(a,1),1,2))
    aconc = ConcatDiskArray(ad)

    cs1 = DiskArrays.eachchunk(aconc)
    @test cs1 isa DiskArrayTools.IrregularGridChunks{3}
    @test cs1.cs[1] == [1:1,2:3,4:6,7:10,11:15,16:21,22:28,29:36,37:45,46:55]
    @test cs1.cs[2] == [i:i for i in 1:15]
    @test cs1.cs[3] == [1:2,3:4]

    @test isapprox(sum(aconc),sum(aar))
#    This will be fixed when broadcast can deal with irregular chunks
    @test_broken copy(aconc .+ ones(1,15)) == aar .+ ones(1,15)
    @test length(ad[3,5].getindex_list) == 10

    @test isapprox.(sum(aconc,dims=1), sum(aar,dims=1)) |> all

    @test length(ad[3,5].getindex_list)==20
end

@testset "Random chunks" begin
    a = [rand(i,j,4) for i = 1:10, j=1:5]
    a1 = [cat(x..., dims = 1) for x in eachcol(a)]
    aar = cat(a1...,dims = 2)
    ad = [_DiskArray(a[i,j], chunksize = (iseven(j) ? 2 : 1,j,4)) for i in 1:10, j in 1:5]
    aconc = ConcatDiskArray(ad)

    cs1 = DiskArrays.eachchunk(aconc)
    allchunks = collect(cs1)

    @test length(allchunks) == 225
    @test isapprox(sum(aconc),sum(aar))
    # This will be fixed when broadcast can deal with irregular chunks
    @test_broken copy(aconc .+ ones(1,15)) == aar .+ ones(1,15)
    @test ad[3,5].getindex_list == [(1:1, 1:5, 1:4), (2:2, 1:5, 1:4),(3:3, 1:5, 1:4)]
    @test ad[3,4].getindex_list == [(1:2, 1:4, 1:4), (3:3, 1:4, 1:4)]

    @test isapprox.(sum(aconc,dims=1), sum(aar,dims=1)) |> all

    @test length(ad[3,5].getindex_list)==6
    @test length(ad[3,4].getindex_list)==4
end

@testset "Regular chunks" begin
    a = [rand(i,j,4) for i = 1:10, j=1:5]
    a1 = [cat(x..., dims = 1) for x in eachcol(a)]
    aar = cat(a1...,dims = 2)
    ad = [_DiskArray(a[i,j], chunksize = (1,1,4)) for i in 1:10, j in 1:5]
    aconc = ConcatDiskArray(ad)

    cs1 = DiskArrays.eachchunk(aconc)

    @test isapprox(sum(aconc),sum(aar))
    @test copy(aconc .+ 1) == aar .+ 1
    @test length(ad[3,5].getindex_list) == 30
    @test length(ad[3,4].getindex_list) == 24

    @test isapprox.(sum(aconc,dims=1), sum(aar,dims=1)) |> all

    @test length(ad[3,5].getindex_list)==45
    @test length(ad[3,4].getindex_list)==36
end




end