"""
    ∩

Compute the intersection between a ray and a geometry. The intersection
can br `nothing` if the ray does not intersect the geometry, a single
point or a tuple of two points.
"""
∩

function ∩(ray::Ray, plane::Plane)
    if ray.direction ⟂ plane.normal
        return ray.origin∈plane ? ray : nothing
    end

    # solve r = origin + α * direction & (r - p₀) * normal = 0
    α = dot(plane.p₀ - ray.origin, plane.normal) / dot(ray.direction, plane.normal)
    if α<0. return nothing end #ray defined only for α≥0
    return ray(α)
end

function ∩(ray::Ray, sphere::Sphere)
    # solve quadratic equation r = origin + α * direction & (r - p₀)^2 = R^2
    a, b, c = dot(ray.direction, ray.direction), 2*dot(ray.direction, ray.origin-sphere.p₀), dot(ray.origin-sphere.p₀,ray.origin-sphere.p₀)-sphere.R^2
    Δ² = b^2 - 4*a*c #squared descriminant
    if Δ² < 0
        return nothing
    elseif abs(Δ²) < 1000*eps(eltype(Δ²))
        α = -b/2a
        return α>=0. ? ray(α) : nothing
    else
        Δ = sqrt(Δ²)
        α₁, α₂ = (-b + Δ)/2a, (-b - Δ)/2a 
        if α₂>=0
            return ray(α₂), ray(α₁)
        elseif α₁>=0
            return ray(α₁)
        else
            return nothing
        end
    end
end

function ∩(ray::Ray, b::Box)
    #solve r = origin + α * direction & vmin ≤ r ≤ vmax 
    α₁, α₂ = (b.vmin - ray.origin) ./ ray.direction, (b.vmax - ray.origin) ./ ray.direction

    #see https://www.cs.cornell.edu/courses/cs4620/2013fa/lectures/03raytracing1.pdf
    α_enter = min.(α₁, α₂)
    α_exit = max.(α₁, α₂)

    α_interval = maximum(α_enter), minimum(α_exit)

    if α_interval[1] > α_interval[2] return nothing
    elseif α_interval[1] < 0.
        if α_interval[2] < 0. return nothing
        else return ray(α_interval[2]) end
    else return ray(α_interval[1]), ray(α_interval[2]) end
end