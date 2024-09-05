"""
    struct Photon

A representation of a photon in a medium with a position `r`, a direction `u` and a weight `weight`.
"""
struct Photon
    "photon position"
    r::SVector{3,Float64}
    "photon direction"
    u::SVector{3,Float64}
    "photon weight"
    weight::Float64

    Photon(r::SVector{3,Float64}, u::SVector{3,Float64}, weight=1.) = new(r, u, weight)
end

Photon(r::Union{AbstractVector,Tuple}, u::Union{AbstractVector,Tuple}, weight=1.) = Photon(SVector(Tuple(r)), SVector(Tuple(u)), weight)

"""
    advance(p::Photon, distance::Float64)

Advance a photon `p` by a distance `distance` in the direction of the photon `p.u`.

# Arguments
- `p::Photon`: the photon to advance
- `distance::Float64`: the distance to advance the photon

# Returns
- a new photon advanced by `distance` in the direction of `p.u`.
"""
@inline function advance(p::Photon, distance::Float64)
    Photon(p.r + distance * p.u, p.u, p.weight)
end

"""
    reflect(p::Photon, normal)

Reflect a photon `p` on a surface with normal `normal` using the law of reflection.

# Arguments
- `p::Photon`: the photon to reflect
- `normal::SVector{3,Float64}`: the normal of the surface

# Returns
- a new photon reflected on the surface
"""
@inline function reflect(p::Photon, normal)
    cosθ1 = -dot(normal, p.u)
    u_reflect = p.u + 2cosθ1*normal
    Photon(p.r, u_reflect, p.weight)
end

@inline function snell_angles(p::Photon, normal, n1::Float64, n2::Float64)
    cosθ1 = -dot(normal, p.u)
    sqrt_body = 1 - (n1/n2)^2 * (1-cosθ1^2)
    if sqrt_body < 0. #pass critical angle
        return cosθ1, nothing
    else 
        cosθ2 = sqrt(sqrt_body)
        return cosθ1, cosθ2
    end
end

"""
    refract(p::Photon, normal, n1::Float64, n2::Float64)

Refract a photon `p` on a surface with normal `normal` where the photon
is coming from a medium with refractive index `n1` to
a medium with refractive index `n2` using Snell's law.

# Arguments
- `p::Photon`: the photon to refract
- `normal::SVector{3,Float64}`: the normal of the surface
- `n1::Float64`: the refractive index of the medium where the photon is coming from
- `n2::Float64`: the refractive index of the medium where the photon is going to

# Returns
- a new photon refracted on the surface
"""
@inline function refract(p::Photon, normal, n1::Float64, n2::Float64)
    cosθ1, cosθ2 = snell_angles(p, normal, n1, n2)
    if isnothing(cosθ2) return nothing end #pass critical angle
    u_refract = (n1/n2)*p.u + ((n1/n2)cosθ1 - cosθ2)*normal
    Photon(p.r, u_refract, p.weight)
end

"""
    unpolarized_reflectance(p::Photon, normal, n1::Float64, n2::Float64)

Compute the unpolarized reflectance (Fresnel reflectance) of a photon `p` on a surface with normal `normal` where the photon
is coming from a medium with refractive index `n1` to a medium with refractive index `n2`.

# Arguments
- `p::Photon`: the photon to reflect
- `normal::SVector{3,Float64}`: the normal of the surface
- `n1::Float64`: the refractive index of the medium where the photon is coming from
- `n2::Float64`: the refractive index of the medium where the photon is going to

# Returns
- the Fresnel reflectance of the photon
"""
@inline function unpolarized_reflectance(p::Photon, normal, n1::Float64, n2::Float64)
    cosθ1, cosθ2 = snell_angles(p, normal, n1, n2)
    if isnothing(cosθ2) return 1. end #pass critical angle = always reflect
    R = 0.5 * ( ((n1*cosθ1-n2*cosθ2)/(n1*cosθ1+n2*cosθ2))^2 + ((n2*cosθ1-n1*cosθ2)/(n2*cosθ1+n1*cosθ2))^2 )
    R
end

"""
    spin(p::Photon, g)

Spin a photon `p` with a Henyey-Greenstein coefficient `g` in a new direction.

# Arguments
- `p::Photon`: the photon to spin
- `g::Float64`: the Henyey-Greenstein coefficient

# Returns
- a new photon with a new random direction
"""
function spin(p::Photon, g)
    if abs(g) < 1000*eps(eltype(g))
        cosθ = 2*rand() - 1
    else
        cosθ = (1 + g^2 - (1-g^2)^2/(1-g+2*g*rand())^2) / 2g
    end
    ψ = 2pi*rand()
    u_rotated = rotate_vector(p.u, cosθ, ψ)
    return Photon(p.r, u_rotated, p.weight)
end

"""
    drop(p::Photon, distnace, μₐ)

Decrease the weight of a photon `p` by a distance `distance` as it is
partially absorbed in the medium with absorption coefficient `μₐ`.

# Arguments
- `p::Photon`: the moving photon
- `distance::Float64`: the distance the photon has traveled
- `μₐ::Float64`: the absorption coefficient of the medium

# Returns
- a new photon with a decreased weight by a factor of `exp(-distance*μₐ)`
"""
@inline function drop(p::Photon, distnace, μₐ)
    return Photon(p.r, p.u, p.weight * exp(-distnace*μₐ))
end

"""
    roulette(p::Photon, chance)

Apply Russian roulette to a photon `p` with a chance `chance` to kill the photon.

# Arguments
- `p::Photon`: the photon to apply Russian roulette
- `chance::Float64`: the chance to kill the photon

# Returns
- a new photon with a increased weight by a factor of `1/chance` if the photon survives, otherwise a photon with weight 0.
"""
@inline roulette(p::Photon, chance) = rand()<chance ? Photon(p.r, p.u, p.weight/chance) : Photon(p.r, p.u, 0.)
