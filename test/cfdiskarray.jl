@testset "CF DiskArray" begin
for T in [Int8, Union{Int8, Missing}]
a = _DiskArray(zeros(T,10,10), chunksize = (2,5))
a2 = CFDiskArray(a,Dict("missing_value"=>-128, "add_offset"=>10, "scale_factor"=>1/128))

@test eltype(a2) == Union{Float64,Missing}

a2[1,:] .= 9.5
a2[2,:] .= 10.5
a2[3,:] .= missing
@test all(isequal(-64), a.parent[1,:])
@test all(isequal(64), a.parent[2,:])
@test all(isequal(-128), a.parent[3,:])
@test all(iszero, a.parent[4:end,:])
@test a2[1,:] == fill(9.5,10)
@test a2[2,:] == fill(10.5,10)
@test all(ismissing,a2[3,:])
end
end
