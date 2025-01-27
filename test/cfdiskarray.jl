@testset "CF DiskArray" begin
a = _DiskArray(zeros(Int8,10,10), chunksize = (2,5))
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

b = _DiskArray([0.0, missing, 2.0, NaN], chunksize=(3,))

b1 = CFDiskArray(b,Dict("missing_value"=>NaN, "add_offset"=>10, "scale_factor"=>2.0))
@test eltype(b1) <: Union{Float64,Missing}
@test all(isequal.(b1[:],[10.0,missing,14.0,missing]))
b1[3:4] = [12.0, 14.0]
b1[1:2] = [missing,missing]
@test isequal.(b1[:],[missing,missing,12.0,14.0]) == [true,true,true,true]
@test isnan(b1.a[2])

b = _DiskArray([0.0, missing, 2.0, NaN], chunksize=(3,))
b2 = CFDiskArray(b,Dict("add_offset"=>0, "scale_factor"=>1.0))
@test b2 isa _DiskArray

b = _DiskArray([0.0, missing, 2.0, NaN], chunksize=(3,))
b3 = CFDiskArray(b,Dict("add_offset"=>0, "scale_factor"=>1.0))
@test eltype(b3) <: Union{Float64,Missing}
@test all(isequal.(b3[:],[0.0, missing, 2.0, NaN]))
b3[3:4] = [12.0, 14.0]
b3[1:2] = [missing,missing]
@test all(isequal.(b1[:],[missing,missing,12.0,14.0]))


b = CFDiskArray(fill(1f0, 3, 3, 3), Dict("add_offset" => 0.0, "scale_factor" => 1.0))
@test eltype(b) == Float32

b = CFDiskArray(fill(1f0, 3, 3, 3),
                Dict("add_offset" => 0.0,
                     "scale_factor" => 1.0,
                     "missing_value" => NaN))
@test eltype(b) == Union{Float32, Missing}

b = CFDiskArray([1f0 ,missing],
                Dict("add_offset" => 0.0,
                     "scale_factor" => 1.0,
                     "missing_value" => NaN))
@test eltype(b) == Union{Float32, Missing}

b = CFDiskArray([1.0f0, missing],
    Dict("add_offset" => 0.0,
        "scale_factor" => 1.0))
@test eltype(b) == Union{Float32,Missing}


b = CFDiskArray([1.0f0, missing], Dict())
@test eltype(b) == Union{Float32,Missing}

b = CFDiskArray([1.0f0, 1.0f0], Dict("missing_value" => NaN))
@test eltype(b) == Union{Float32,Missing}

b = CFDiskArray([1.0f0, 2.0f0],
    Dict("add_offset" => 0.0,
        "scale_factor" => 1.0))
@test eltype(b) == Float32

b = CFDiskArray([10, -9999],
    Dict("scale_factor" => Float16(0.1),
        "missing_value" => -9999))
@test eltype(b) == Union{Missing, Float16}
@test ismissing(b[2])

# CF conventions prescribe a "_FillValue field"
b = CFDiskArray([1.0f0, 2.0f0],
                Dict("_FillValue" => NaN32))
@test eltype(b) == Float32
end
