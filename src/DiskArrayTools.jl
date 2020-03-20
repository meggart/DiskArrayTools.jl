module DiskArrayTools
import DiskArrays: AbstractDiskArray, eachchunk, haschunks, Chunked, estimate_chunksize, GridChunks, findints, readblock!, writeblock!

export DiskArraySTack, diskstack
struct DiskArrayStack{T,N,M,NO}<:AbstractDiskArray{T,N}
    arrays::Array{M,NO}
end
function diskstack(a::Array{M,N}) where {N,M<:AbstractArray{T,NO}} where {T,NO}
    length(unique(size.(a))) == 1 || error("All arrays in the stack must have the same size")
    DiskArrayStack{T,N+NO,M,N}(a)
end
Base.size(r::DiskArrayStack) = (size(r.arrays[1])...,size(r.arrays)...)
haschunks(a::DiskArrayStack) = haschunks(a.arrays[1])
iscompressed(a::DiskArrayStack) = any(iscompressed,a.arrays)
function eachchunk(a::DiskArrayStack{<:Any,<:Any,<:Any,NO}) where NO
    iterold = eachchunk(a.arrays[1])
    cs = (iterold.chunksize...,ntuple(one,NO)...)
    co = (iterold.offset...,ntuple(zero,NO)...)
    DiskArrays.GridChunks(a,cs,offset=co)
end

function DiskArrays.readblock!(a::DiskArrayStack{<:Any,N,<:Any,NO},aout,i::AbstractVector...) where {N,NO}

  innerinds = i[1:(N-NO)]

  outerinds = i[(N-NO+1):N]
  innercolon = map(_->(:), innerinds)
  iiter = CartesianIndices(outerinds)
  inum  = CartesianIndices(size(iiter))
  foreach(zip(iiter,inum)) do (iouter,iret)
    arnow = a.arrays[iouter]
    if isa(arnow, AbstractDiskArray)
      readblock!(a.arrays[iouter],view(aout,innercolon...,iret.I...),innerinds...)
    else
      aout[innercolon...,iret.I...] = arnow[innerinds...]
    end
  end
  nothing
end

function DiskArrays.writeblock!(a::DiskArrayStack{<:Any,N,<:Any,NO},v,i::AbstractVector...) where {N,NO}
  innerinds = i[1:(N-NO)]

  outerinds = i[(N-NO+1):N]
  innercolon = map(_->(:), innerinds)
  iiter = CartesianIndices(outerinds)
  inum  = CartesianIndices(size(iiter))
  foreach(zip(iiter,inum)) do (iouter,iret)
    arnow = a.arrays[iouter]
    if isa(arnow, AbstractDiskArray)
      writeblock!(a.arrays[iouter],view(v,innercolon...,iret.I...),innerinds...)
    else
      arnow[innerinds...] = aout[innercolon...,iret.I...]
    end
  end
  nothing
end

function Base.view(a::DiskArrayStack{<:Any,N,<:Any,NO},i...) where {N,NO}
    iinner = i[1:(N-NO)]
    ashort = map(view(a.arrays,i[(N-NO+1):N]...)) do ai
        view(ai,iinner...)
    end
    if ndims(ashort)==0
      ashort[]
    else
      diskstack(ashort)
    end
end

using DiskArrays: AbstractDiskArray, writeblock!, toRanges, DiskArrays, GridChunks
import DiskArrays: readblock!, haschunks, eachchunk
using Interpolations
using Base.Iterators: product
using OffsetArrays: OffsetArray
struct InterpolatedDiskArray{T,N,A<:AbstractArray{T,N},I,O,BC,CS<:Union{Nothing,NTuple{N,Int}}} <: AbstractDiskArray{T,N}
    a::A
    newinds::I
    meth::O
    bc::BC
    chunksize::CS
end
allmeths(order::Tuple,newinds) = map(order,newinds) do o,ni
    ni === nothing ? NoInterp() : BSpline(o)
end
function allmeths(order, newinds)
    allorders = map(newinds) do _
        order
    end
    allmeths(allorders,newinds)
end
function InterpolatedDiskArray(a::AbstractArray,chunksize,newinds...; order=Linear(), bc=Flat())
    eltype(a) <: Union{Real,Complex} || error("Can only interpolate real or complex values")
    s = size(a)
    ni2 = map((s,i)->i===nothing ? (1:s) : i,s,newinds)
    me = allmeths(order,newinds)
    InterpolatedDiskArray(a,ni2,me,bc,chunksize)
end
round_readinds(s,r) = max(1,floor(Int,first(r))):min(s,ceil(Int,last(r)))
Base.size(a::InterpolatedDiskArray) = map(length,a.newinds)
haschunks(a::InterpolatedDiskArray) = a.chunksize===nothing
eachchunk(a::InterpolatedDiskArray) = GridChunks(a,a.chunksize)



function fremap(xout,xin,newinds,ipmeth;bc = Flat())
  interp = extrapolate(interpolate(xin, ipmeth), bc)
  for (i,ic) in enumerate(Iterators.product(newinds...))
    xout[i]=interp(ic...)
  end
end
function readblock!(a::InterpolatedDiskArray, aout, i::AbstractUnitRange...)
    parentranges = map((ir,ip)->ip[ir],i,a.newinds)
    rr = round_readinds.(size(a.a),parentranges)
    atemp = a.a[rr...]
    meth = map((s,m)->s==1 ? NoInterp() : m,size(atemp),a.meth)
    atemp2 = OffsetArray(atemp,map(r->(first(r)-1),rr))
    if all(isequal(NoInterp()),meth)
        aout .= atemp
    else
        fremap(aout,atemp2,parentranges,rr,meth,bc=a.bc)
    end
end


end # module
