"""
 struct Scene

A collection of solids in the scene. The order of the solids is important,
as the earliest solids "hide" the later ones.
"""    
struct Scene
    "The ordered collections of solids in the scene, ordered by priority"
    solids::Vector{Solid}
    Scene() = new([Outside()])
end

∈(p::Photon, s::Scene) = any(p.r∈sol.geometry for sol in s.solids)

"""
    add_solid!(s::Scene, solid::Solid)

Add a solid to the scene.
"""
add_solid!(s::Scene, solid::Solid) = pushfirst!(s.solids, solid)

"""
    add_solid!(s::Scene, geometry, material)

Add a solid to the scene with the given geometry and material.
"""
add_solid!(s::Scene, geometry, material) = add_solid!(s, Solid(geometry, material))

"""
    Scene(solid_collection...)

Create a scene from a collection of solids.
"""
function Scene(solid_collection)
    scene = Scene()
    pushfirst!(scene.solids, solid_collection...)
    scene
end

"""
    find_container(p::Photon, s::Scene)

Find the solid that contains the photon. The solid is the first solid in the
scene's list that contains the photon.
"""
find_container(p::Photon, s::Scene) = first(sol for sol in s.solids if p.r∈sol.geometry)