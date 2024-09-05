"""
    struct Sphere{D,T} <: Geometry{D,T}

A sphere geometry object in the scene.
"""
struct Sphere{D,T} <: Geometry{D,T}
    "The center of the sphere"
    p₀::SVector{D,T}
    "The radius of the sphere"
    R::T
end

Sphere(p₀::NTuple{D,T}, R::T) where {D,T} = Sphere{D,T}(p₀, R)
Sphere(p₀::Union{AbstractVector,Tuple}, R) = Sphere(Tuple(p₀), R)

∈(x::AbstractVector, sphere::Sphere) = dot(x-sphere.p₀,x-sphere.p₀) < sphere.R^2 || isapprox(dot(x-sphere.p₀,x-sphere.p₀), sphere.R^2; atol=100*eps(Float64))
