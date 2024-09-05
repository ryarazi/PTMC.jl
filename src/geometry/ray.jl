"""
    struct Ray{D,T} <: Geometry{D,T}

A ray geometry object in the scene.
"""
struct Ray{D,T} <: Geometry{D,T}
    "The origin of the ray"
    origin::SVector{D,T}
    "The direction vector of the ray"
    direction::SVector{D,T}
end

Ray(origin::NTuple{D,T}, direction::NTuple{D,T}) where {D,T} = Ray{D,T}(origin, direction)
Ray(origin::Union{AbstractVector,Tuple}, direction::Union{AbstractVector,Tuple}) = Ray(Tuple(origin), Tuple(direction))

@inline function (r::Ray)(α)
    # @assert α>=0.
    r.origin + α*r.direction
end