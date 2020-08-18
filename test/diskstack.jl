@testset "DiskStack" begin
ars = [rand(10,5) for i = 1:3]
cars = cat(ars...,dims=3)
dars = map(a->_DiskArray(a,chunksize=(5,2)),ars)
ds = diskstack(dars)
@test size(ds) == (10,5,3)
@test ds isa DiskArrayStack{Float64,3,_DiskArray{Float64,2,Array{Float64,2}},1}
@test isapprox(sum(ds), sum(cars))
@test all(i->getindex_count(i) == 6,dars)
@test ds[2:3,3:4,:] == cars[2:3,3:4,:]
@test all(i->getindex_count(i) == 7,dars)
v = view(ds,:,:,2)
@test v isa DiskArrays.SubDiskArray{Float64,2}
@test v[:,:] == ars[2]
@test getindex_count(dars[2]) == 8
v2 = view(ds,:,:,[3,1])
@test v2[:,:,:] == cars[:,:,[3,1]][:,:,:]
@test all(i->getindex_count(i) == 8,dars)
end
