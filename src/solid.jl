"""
    Solid

A solid is a geometry with a material. It is the basic building block of a scene.
"""
struct Solid
    geometry::Geometry{3,Float64}
    material::Turbid
    name::String
    Solid(geometry, material; name="") = new(geometry, material, name)
end

"""
    Outside()

Return a solid that represents the outside of the scene.
"""
Outside() = Solid(Box((-Inf,-Inf,-Inf), (Inf,Inf,Inf)), Turbid(), name="Outside")

∈(p::Photon, s::Geometry{3,Float64}) = p.r∈s
∈(p::Photon, s::Solid) = p∈s.geometry

"""
    exit_solid(p::Photon, s::Solid)

Return the point where the photon exits the solid and the normal to the surface at that point.

# Arguments
- `p::Photon`: the photon moving through the solid
- `s::Solid`: the solid the photon is moving through

# Returns
- `(boundary_intersection, normal)`: the point where the photon exits the solid and the normal to the surface at that point
"""
exit_solid(p::Photon, s::Solid) = exit_solid(p, s.geometry)

function exit_solid(p::Photon, b::Box{3,Float64})
    # @assert p∈b
    ray = Ray(p.r, p.u)
    boundary_intersection::SVector{3, Float64} = ray∩b
    normal = (abs.(boundary_intersection-b.vmin) .< 1e5*eps(Float64)) - (abs.(boundary_intersection-b.vmax) .< 1e5*eps(Float64))
    # @assert sum(abs.(normal))==1
    return boundary_intersection, normal
end