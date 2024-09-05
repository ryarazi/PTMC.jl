"""
    struct Box{D,T} <: Geometry{D,T}

A box geometry object in the scene.
"""
struct Box{D,T} <: Geometry{D,T}
    vmin::SVector{D,T}
    vmax::SVector{D,T}

    function Box{D,T}(vmin, vmax) where {D,T}
        @assert all(vmin .<= vmax)
        new(vmin, vmax)
    end
end

Box(vmin::NTuple{D,T}, vmax::NTuple{D,T}) where {D,T} = Box{D,T}(vmin, vmax)
Box(vmin::Union{AbstractVector,Tuple}, vmax::Union{AbstractVector,Tuple}) = Box(Tuple(vmin), Tuple(vmax))

#overload Base methods
minimum(b::Box) = b.vmin
maximum(b::Box) = b.vmax
extrema(b::Box) = b.vmin, b.vmax
∈(x::AbstractVector, b::Box) = all(@. ((b.vmin < x) || isapprox(b.vmin, x; atol=100*eps(Float64))) & ((x < b.vmax) || isapprox(b.vmax, x; atol=100*eps(Float64))))