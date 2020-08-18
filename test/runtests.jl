using DiskArrayTools
using Test
using DiskArrays: DiskArrays, ReshapedDiskArray, PermutedDiskArray, AbstractDiskArray

#Define a data structure that can be used for testing
struct _DiskArray{T,N,A<:AbstractArray{T,N}} <: AbstractDiskArray{T,N}
  getindex_list::Vector{NTuple{N,AbstractUnitRange}}
  setindex_list::Vector{NTuple{N,AbstractUnitRange}}
  parent::A
  chunksize::NTuple{N,Int}
end
_DiskArray(a;chunksize=size(a)) = _DiskArray(NTuple{ndims(a),AbstractUnitRange}[],NTuple{ndims(a),AbstractUnitRange}[],a,chunksize)
Base.size(a::_DiskArray) = size(a.parent)
DiskArrays.haschunks(::_DiskArray) = DiskArrays.Chunked()
DiskArrays.eachchunk(a::_DiskArray) = DiskArrays.GridChunks(a,a.chunksize)
getindex_count(a::_DiskArray) = length(a.getindex_list)
setindex_count(a::_DiskArray) = length(a.setindex_list)
trueparent(a::_DiskArray) = a.parent
getindex_count(a::ReshapedDiskArray) = getindex_count(a.parent)
setindex_count(a::ReshapedDiskArray) = setindex_count(a.parent)
trueparent(a::ReshapedDiskArray) = trueparent(a.parent)
getindex_count(a::PermutedDiskArray) = getindex_count(a.a.parent)
setindex_count(a::PermutedDiskArray) = setindex_count(a.a.parent)
trueparent(a::PermutedDiskArray{T,N,<:PermutedDimsArray{T,N,perm,iperm}}) where {T,N,perm,iperm} = permutedims(trueparent(a.a.parent),perm)
function DiskArrays.readblock!(a::_DiskArray,aout,i::AbstractUnitRange...)
  ndims(a) == length(i) || error("Number of indices is not correct")
  all(r->isa(r,AbstractUnitRange),i) || error("Not all indices are unit ranges")
  #println("reading from indices ", join(string.(i)," "))
  push!(a.getindex_list,i)
  aout .= a.parent[i...]
end
function DiskArrays.writeblock!(a::_DiskArray,v,i::AbstractUnitRange...)
  ndims(a) == length(i) || error("Number of indices is not correct")
  all(r->isa(r,AbstractUnitRange),i) || error("Not all indices are unit ranges")
  #println("Writing to indices ", join(string.(i)," "))
  push!(a.setindex_list,i)
  view(a.parent,i...) .= v
end

include("diskstack.jl")
include("interpolate.jl")
include("cfdiskarray.jl")
