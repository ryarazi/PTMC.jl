"""
    struct Plane{D,T} <: Geometry{D,T}

A plane geometry object in the scene.
"""
struct Plane{D,T} <: Geometry{D,T}
    "A point on the plane"
    p₀::SVector{D,T}
    "The normal to the plane"
    normal::SVector{D,T}
end

Plane(p₀::NTuple{D,T}, normal::NTuple{D,T}) where {D,T} = Plane{D,T}(p₀, normal)

∈(x::AbstractVector, plane::Plane) = (x⋅plane.normal ≈ plane.p₀⋅plane.normal)