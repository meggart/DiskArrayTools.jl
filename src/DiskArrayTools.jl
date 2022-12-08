module DiskArrayTools
import DiskArrays: AbstractDiskArray, eachchunk, haschunks, Chunked,
estimate_chunksize, GridChunks, findints, readblock!, writeblock!, 
RegularChunks, IrregularChunks, ChunkType, approx_chunksize, chunktype_from_chunksizes
using Interpolations
using IterTools: imap
using Base.Iterators: product
using OffsetArrays: OffsetArray


export DiskArrayStack, diskstack, ConcatDiskArray, CFDiskArray
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
function eachchunk(a::DiskArrayStack)
    oldchunks = eachchunk(a.arrays[1]).chunks
    newchunks = map(s->RegularChunks(1,0,s),size(a.arrays))
    GridChunks(oldchunks..., newchunks...)
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

function readblock!(a::ResampledDiskArray, aout, i::AbstractUnitRange...)
  parentranges = map((ir,ip)->ip[ir],i,a.newinds)
  rr = get_readinds(a,parentranges)
  atemp = a.a[rr...]
  atemp2 = OffsetArray(atemp,map(r->(first(r)-1),rr))
  resample_disk(a,aout,atemp2,parentranges)
end



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
  allnearest = all(i->isa(i,Union{Nothing,BSpline{<:Constant}}),a.meth)
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

function aggregate_chunks(cs,agg)
  lcs = cumsum(length.(cs))
  lagg = cumsum(length.(agg))
  inds = [searchsortedfirst(lagg,lc) for lc in lcs]
  chunktype_from_chunksizes(diff([0;inds]))
end
interpret_aggsize(_, x::ChunkType) = x
interpret_aggsize(sparent,x) = RegularChunks(x,0,sparent)


struct AggregatedDiskArray{T,N,A<:AbstractArray{T,N},I,F,CS<:Union{Nothing,GridChunks{N}}} <: ResampledDiskArray{T,N}
  a::A
  newinds::I
  aggfun::F
  chunksize::CS
  dimkw::Bool
end
Base.size(a::AggregatedDiskArray) = map(length,a.newinds)
eachchunk(a::AggregatedDiskArray) = a.chunksize
function AggregatedDiskArray(parent,aggsizes,aggfun;chunksizes=nothing,dimkw=true)
  newinds = GridChunks(interpret_aggsize.(size(parent),aggsizes)...)
  if chunksizes === nothing
      chunksizes = GridChunks(aggregate_chunks.(eachchunk(parent).chunks,newinds.chunks)...)
  end
  AggregatedDiskArray(parent,newinds.chunks,aggfun,chunksizes,dimkw)
end


get_readinds(::AggregatedDiskArray, parentranges) = map(r->first(first(r)):last(last(r)),parentranges)
function resample_disk(a::AggregatedDiskArray,aout,atemp2,parentranges)
  if a.dimkw
    resample_disk_dimkw(a,aout,atemp2,parentranges)
  else
    resample_disk_loop(a,aout,atemp2,parentranges)
  end
end

#Generic fallback
function resample_disk_loop(a::AggregatedDiskArray,aout,atemp2,parentranges)
  allr = collect(Iterators.product(parentranges...))
  for ic in eachindex(allr)
      ii = allr[ic...]
      aout[ic...] = a.aggfun(atemp2[ii...])
  end
end

function resample_disk_dimkw(a::AggregatedDiskArray,aout,atemp2,parentranges)
  allpure = map(parentranges) do pr
    all(i->length(i)==1,pr) ? [Colon()] : pr
  end
  alloutinds = map(allpure) do pr
    isa(pr[1],Colon) ? pr : 1:length(pr)
  end
  dimkw = (findall(i->i[1] !== Colon(), allpure)...,)
  allr = Iterators.product(allpure...)
  outr = Iterators.product(alloutinds...)
  for (ii,iout) in zip(allr, outr)
      aout[iout...] = a.aggfun(atemp2[ii...],dims=dimkw)
  end
end


remove_missing(::Type{T}) where T <: Union{Missing, T2} where T2 = T2

#Use of Sentinel missing value
struct CFDiskArray{T,N,MT,P,OT} <: AbstractDiskArray{T,N}
    a::P
    mv::MT
    add_offset::OT
    scale_factor::OT
end
function CFDiskArray(a::AbstractArray{T}, attr::Dict) where T
  mv = get(attr,"missing_value", nothing)
  offs, sc = if haskey(attr, "add_offset") || haskey(attr, "scale_factor")
    _offs = get(attr, "add_offset", false)
    _sc = get(attr, "scale_factor", true)
    T2 = remove_missing(T)
    if _offs isa AbstractFloat && _sc isa AbstractFloat && T2 <: AbstractFloat
      T2(_offs), T2(_sc)
    else
      promote(_offs, _sc)
    end
  else
    zero(T), one(T)
  end
  #Short track if there is nothing to do
  if mv === nothing && offs==zero(offs) && sc == one(sc)
    return a
  end
  # Now check if element type of a is ok to use
  T_pure = typeof(offs)
  if T >: Missing 
    if mv !== nothing
      @warn "Trying to construct a CFDiskArray from an array that already contains missings. In case this comes from a Zarr dataset, consider opening the dataset with `fill_as_missing=false`"
    else
      T_pure = Union{T_pure,Missing}
    end
  end
  S,mv = if mv === nothing
    T_pure,mv
  else
    Union{T_pure,Missing},convert(T_pure,mv)
  end
  CFDiskArray{S,ndims(a),typeof(mv),typeof(a),typeof(offs)}(a, mv, offs, sc)
end

Base.size(a::CFDiskArray, args...) = size(a.a, args...)
haschunks(a::CFDiskArray) = haschunks(a.a)
eachchunk(a::CFDiskArray) = eachchunk(a.a)
iscompressed(a::CFDiskArray) = iscompressed(a.a)

function readblock!(a::CFDiskArray, aout, r::AbstractVector...)
    mv = a.mv
    sc,offs = a.scale_factor, a.add_offset
    broadcast!(scaleoffs, aout, a.a[r...],mv,sc,offs)
    nothing
end
checkmiss(x,mv) = isequal(x,mv)
checkmiss(_,::Nothing) = false
scaleoffs(x,mv,sc,offs) = checkmiss(x,mv) ? missing : x*sc+offs
function writeblock!(a::CFDiskArray, v, r::AbstractVector...)
    mv = a.mv
    sc,offs = a.scale_factor, a.add_offset
    a.a[r...] = broadcast(scaleoffsinv, v,mv,sc,offs)
    nothing
end
scaleoffsinv(x,mv::Integer,sc,offs) = ismissing(x) ? mv : round(typeof(mv),(x-offs)/sc)
scaleoffsinv(x,mv,sc,offs) = ismissing(x) ? mv : ((x-offs)/sc)
scaleoffsinv(x,::Nothing,sc,offs) = (x-offs)/sc

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
haschunks(::ConcatDiskArray) = Chunked()
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
        aout[outer_range...] = a.parents[cI][array_range...]
    end
end
function writeblock!(a::ConcatDiskArray, aout, inds::AbstractUnitRange...)
  error("No method yet for writing into a ConcatDiskArray")
end

function mergechunks(a::RegularChunks, b::RegularChunks)
  if a.s==0 || (a.cs == b.cs && length(last(a))==a.cs)
    RegularChunks(a.cs, a.offset,a.s+b.s)
  else
    mergechunks_irregular(a,b)
  end
end

mergechunks(a::ChunkType, b::ChunkType) = mergechunks_irregular(a,b)
function mergechunks_irregular(a, b)
  IrregularChunks(chunksizes = filter(!iszero,[length.(a); length.(b)]))
end

#This function will work for DiskArrayStack as well, we only need a fallback_chunksize function
function eachchunk(aconc::ConcatDiskArray)
  N = ndims(aconc)
  s = size(aconc)
  oldchunks = eachchunk.(aconc.parents)
  newchunks = ntuple(N) do i
    sliceinds = Base.setindex(map(_ -> 1,s),:,i)
    v = map(c->c.chunks[i],oldchunks[sliceinds...])
    init = RegularChunks(approx_chunksize(first(v)),0,0)
    reduce(mergechunks, v, init=init)
  end
  GridChunks(newchunks...)
end

end # module
