
"""
    simulate_particle(source::Source, scene::Scene, aggregators = nothing; roulette_threshold = 0.001, roulette_chance = 0.1, δs = 1.e-8)

Simulate a single particle through the scene. The particle is created by sampling the source and then
propagated through the scene until it exits the system or its weight is zero. The particle is then returned.

# Arguments
- `source::Source`: The source from which the particle is sampled.
- `scene::Scene`: The scene through which the particle is propagated.
- `aggregators`: An Aggregator subtype or a collection of Aggregators which collect data during the simulation.
- `roulette_threshold::Float64`: The threshold at which the particle weight is considered small enough to be terminated.
- `roulette_chance::Float64`: The probability that a particle is terminated when its weight is below the threshold.
- `δs::Float64`: The distance to move the particle past the boundary before it is reflected or refracted.

# Returns
- `Photon`: The particle after it has been propagated through the scene.
"""
function simulate_particle(source::Source, scene::Scene, aggregators = nothing; roulette_threshold = 0.001, roulette_chance = 0.1, δs = 1.e-8)  
    p = sample_source(source)
    solid = find_container(p, scene)

    if isa(aggregators, Aggregator) aggregators = [aggregators] end
    if !isnothing(aggregators) for agg in aggregators create_aggregation(agg, p) end end
    
    while p.weight != 0. && solid != Outside()
        s = randexp() # HOP

        #check boundaries and move between solids
        while s != 0.
            pstart = p
            if solid == Outside() break end #photon is out of the systems
            
            boundary_intersection, normal = exit_solid(p, solid)
            d_boundary = sqrt.(sum((p.r - boundary_intersection).^2)) #calculate norm, faster then norm function
            s_boundary = d_boundary * solid.material.μₛ

            if s_boundary > s
                d = s/solid.material.μₛ
                p = advance(p, d)
                p = spin(p, solid.material.g) # SPIN
                p = drop(p, d, solid.material.μₐ) # DROP    
            else
                p = Photon(boundary_intersection, p.u, p.weight) #move photon to the surface boundary
                solid_next = find_container(advance(p, δs), scene)
                isreflect = rand() < unpolarized_reflectance(p, normal, solid.material.n, solid_next.material.n) #decide if reflect or refract
                p = isreflect ? reflect(p, normal) : refract(p, normal, solid.material.n, solid_next.material.n)
                p = advance(p, δs) #make sure the particle is just past the boundary
                p = drop(p, d_boundary, solid.material.μₐ) # DROP
                solid = isreflect ? solid : solid_next #decide on next solid
                # @assert p∈solid.geometry
            end
            s = max(s - s_boundary, 0.)

            if !isnothing(aggregators) for agg in aggregators update_aggregation(agg, pstart, p) end end

            #Roulette Method (see https://fa.bianp.net/blog/2022/russian-roulette/)
            if p.weight < roulette_threshold
                p = roulette(p, roulette_chance)
            end
        end
    end
    
    return p
end