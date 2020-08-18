@testset "Interpolation" begin
using DiskArrayTools: InterpolatedDiskArray
using DiskArrayTools.Interpolations
using DiskArrays: GridChunks

a = reshape(1.0:40,8,5)
ad = _DiskArray(a, chunksize=(4,3))
adi = InterpolatedDiskArray(ad,GridChunks((17,5),(8,3)),0.5:0.5:8.5,nothing)
@test size(adi) == (17,5)
@test adi[15:17,4:5] == [31.5 39.5; 32.0 40.0; 32.0 40.0]
@test adi[1,:] == a[1,:]
@test adi[1:3,1] == [1.0,1.0,1.5]
@test getindex_count(ad) == 3
adi2 = InterpolatedDiskArray(ad,GridChunks((8,2),(4,4)),[1.1,1.2,1.3,1.4,1.6,1.7,1.8,1.9],[1.25,1.75], order = (Constant(),Linear()))
@test adi2[:,:] == [fill(3.0,4) fill(7.0,4); fill(4.0,4) fill(8.0,4)]
ad = _DiskArray(a, chunksize=(4,3))
adi3 = InterpolatedDiskArray(ad,GridChunks((3,3),(2,2)),[1.5,2.5,3.5],[1.25,1.5,1.75], order = (Constant(),Constant()))
@test all(adi3[1:2,1:2] .== a[2,1])
end
