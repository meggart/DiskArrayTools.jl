using DiskArrays, DiskArrayTools
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