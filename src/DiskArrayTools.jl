module DiskArrayTools
import DiskArrays: AbstractDiskArray, eachchunk, haschunks, Chunked,
estimate_chunksize, GridChunks, findints, readblock!, writeblock!
using Interpolations
using Base.Iterators: product
using OffsetArrays: OffsetArray


export DiskArrayStack, diskstack, ConcatDiskArray
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
    GridChunks(a,cs,offset=co)
end

function readblock!(a::DiskArrayStack{<:Any,N,<:Any,NO},aout,i::AbstractVector...) where {N,NO}

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

function writeblock!(a::DiskArrayStack{<:Any,N,<:Any,NO},v,i::AbstractVector...) where {N,NO}
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
      arnow[innerinds...] = v[innercolon...,iret.I...]
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

abstract type ResampledDiskArray{T,N} <: AbstractDiskArray{T,N} end
get_readinds(a::ResampledDiskArray, r) = error("get_readinds not implemented for $(typeof(a))")



struct InterpolatedDiskArray{T,N,A<:AbstractArray{T,N},I,O,BC,CS<:Union{Nothing,GridChunks{N}}} <: ResampledDiskArray{T,N}
    a::A
    newinds::I
    meth::O
    bc::BC
    chunksize::CS
end
fixin(i,s) = max(1,min(s,i))
function round_readinds(s,r,an)
  mi,ma = extrema(r)
  if mi == ma
    m = fixin(round(Int,mi),s)
    m:m
  else
    a = if an
      (round(Int,mi,RoundNearestTiesUp), -round(Int,-ma,RoundNearestTiesUp))
    else
      (floor(Int,mi),ceil(Int,ma))
    end
    fixin(a[1],s):fixin(a[2],s)
  end
end

function get_readinds(a::InterpolatedDiskArray,r)
  allnearest = all(i->isa(i,Union{Nothing,BSpline{Constant}}),a.meth)
  round_readinds.(size(a.a),r,allnearest)
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
    eltype(a) <: Union{Missing,Real,Complex} || error("Can only interpolate real or complex values")
    s = size(a)
    ni2 = map((s,i)->i===nothing ? (1:s) : i,s,newinds)
    me = allmeths(order,newinds)
    InterpolatedDiskArray(a,ni2,me,bc,chunksize)
end
Base.size(a::InterpolatedDiskArray) = map(length,a.newinds)
haschunks(a::InterpolatedDiskArray{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:GridChunks}) = true
haschunks(a::InterpolatedDiskArray{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,Nothing}) = false
eachchunk(a::InterpolatedDiskArray{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:GridChunks}) = a.chunksize

function resample_disk(a::InterpolatedDiskArray,aout,atemp,parentranges)
  meth = map((s,m)->s==1 ? NoInterp() : m,size(atemp),a.meth)
  if all(isequal(NoInterp()),meth)
      aout .= atemp.parent
  else
      fremap(aout,atemp,parentranges,meth,bc=a.bc)
  end
end

function fremap(xout,xin,newinds,ipmeth;bc = Flat())
  interp = extrapolate(interpolate(xin, ipmeth), bc)
  for (i,ic) in enumerate(Iterators.product(newinds...))
    xout[i]=interp(ic...)
  end
end
function readblock!(a::ResampledDiskArray, aout, i::AbstractUnitRange...)
    parentranges = map((ir,ip)->ip[ir],i,a.newinds)
    rr = get_readinds(a,parentranges)
    atemp = a.a[rr...]
    atemp2 = OffsetArray(atemp,map(r->(first(r)-1),rr))
    resample_disk(a,aout,atemp2,parentranges)
end


#Use of Sentinel missing value
struct CFDiskArray{T,N,ST,P<:AbstractArray{ST,N}} <: AbstractDiskArray{Union{T,Missing},N}
    a::P
    mv::ST
    add_offset::T
    scale_factor::T
end
function CFDiskArray(a::AbstractArray{T}, attr::Dict) where T
  mv = get(attr,"missing_value",get(attr,"_FillValue",typemax(T)))
  offs,sc = if haskey(attr,"add_offset") || haskey(attr,"scale_factor")
    offs = get(attr,"add_offset",zero(Float16))
    sc = get(attr,"scale_factor",one(Float16))
    promote(offs,sc)
  else
    zero(T), one(T)
  end
  CFDiskArray(a, mv, offs, sc)
end

Base.size(a::CFDiskArray, args...) = size(a.a, args...)
haschunks(a::CFDiskArray) = haschunks(a.a)
eachchunk(a::CFDiskArray) = eachchunk(a.a)
iscompressed(a::CFDiskArray) = iscompressed(a.a)

function readblock!(a::CFDiskArray, aout, r::AbstractVector...)
    mv = a.mv
    sc,offs = a.scale_factor, a.add_offset
    map!(j->scaleoffs(j,mv,sc,offs), aout, a.a)
    nothing
end
scaleoffs(x,mv,sc,offs) = x==mv ? missing : x*sc+offs
function writeblock!(a::CFDiskArray, v, r::AbstractVector...)
    mv = a.mv
    sc,offs = a.scale_factor, a.add_offset
    map!(j->scaleoffsinv(j,mv,sc,offs), a.a, v)
    nothing
end
scaleoffsinv(x,mv::Integer,sc,offs) = ismissing(x) ? mv : round(typeof(mv),(x-offs)/sc)
scaleoffsinv(x,mv,sc,offs) = ismissing(x) ? mv : ((x-offs)/sc)


struct ConcatDiskArray{T,N,P} <: AbstractDiskArray{T,N}
    parents::P
    startinds::NTuple{N,Vector{Int}}
    size::NTuple{N,Int}
end

function ConcatDiskArray(arrays::AbstractArray)
    #First do some consistency checks
    arrinfo = map(arrays) do a
        eltype(a), size(a)
    end
    et = first(arrinfo)[1]
    nd = length(first(arrinfo)[2])
    sa = size.(arrays)
    all(i->i[1]==et,arrinfo) || error("Arrays don't have the same element type")
    all(i->length(i[2])==nd, arrinfo) || error("Arrays don't have the same dimensions")
    if ndims(arrays)<nd
        arrays = reshape(arrays,size(arrays)...,ones(Int,nd-ndims(arrays))...)
    end
    function othersize(x,id)
        (x[1:id-1]..., x[id+1:end]...)
    end
    si = map(1:nd) do id
        otherdims = ((1:id-1)...,(id+1:nd)...)
        a = size.(arrays)
        a = reduce(a, dims=id, init=ntuple(zero,nd)) do i,j
                if all(iszero,i)
                    j
                elseif othersize(i,id)==othersize(j,id)
                    j
                else
                    error("Dimensions don't match")
                end
        end
        I = fill!(Array{Union{Int,Colon},1}(undef,nd),1)
        I[id] = Colon()
        ar = sa[I...]
        ari = map(i->i[id],ar)
        sl = sum(ari)
        r = cumsum(ari)
        pop!(pushfirst!(r,0))
        r.+1, sl
    end
    ConcatDiskArray{et,nd,typeof(arrays)}(arrays, (getindex.(si,1)...,),(getindex.(si,2)...,))
end

Base.size(a::ConcatDiskArray) = a.size

import DiskArrays: readblock!, writeblock!
using OffsetArrays: OffsetArray
function readblock!(a::ConcatDiskArray, aout, inds::AbstractUnitRange...)
    #Find affected blocks and indices in blocks
    blockinds  = map(inds, a.startinds, size(a.parents)) do i, si, s
        bi1 = max(searchsortedlast(si,first(i)),1)
        bi2 = min(searchsortedfirst(si,last(i)+1)-1,s)
        bi1:bi2
    end
    map(CartesianIndices(blockinds)) do cI
        myar = a.parents[cI]
        mysize = size(myar)
        array_range = map(cI.I, a.startinds, mysize, inds) do ii, si, ms, indstoread
            max(first(indstoread)-si[ii]+1,1) : min(last(indstoread)-si[ii]+1, ms)
        end
        outer_range = map(cI.I, a.startinds, array_range, inds) do ii, si, ar, indstoread
            (first(ar)+si[ii]-first(indstoread)):(last(ar)+si[ii]-first(indstoread))
        end
        aout[outer_range...] .= a.parents[cI][array_range...]
    end
end

end # module
