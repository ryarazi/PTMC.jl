"""
    abstract type Source

An abstract type representing a generic source of particles.
"""
abstract type Source end

"""
    struct Laser <: Source

A type representing a laser source.

# Fields
- `σ::Float64`: The standard deviation of the laser beam.
- `direction::SVector{3,Float64}`: The direction vector of the laser beam.
- `center::SVector{3,Float64}`: The center position of the laser beam.

# Constructor
- `Laser(σ, direction::SVector{3,Float64}, center::SVector{3,Float64})`: Creates a new `Laser` instance with the given standard deviation, direction, and center. The direction vector is normalized.
"""
struct Laser <: Source
    σ::Float64
    direction::SVector{3,Float64}
    center::SVector{3,Float64}

    Laser(σ, direction::SVector{3,Float64}, center::SVector{3,Float64}) = new(σ, normalize(direction), center)
end

"""
    Laser(σ, origin::Union{AbstractVector,Tuple}, direction::Union{AbstractVector,Tuple}) -> Laser

Alternative constructor for `Laser` that accepts the origin and direction as vectors or tuples.

# Arguments
- `σ::Float64`: The standard deviation of the laser beam.
- `origin::Union{AbstractVector,Tuple}`: The origin position of the laser beam.
- `direction::Union{AbstractVector,Tuple}`: The direction vector of the laser beam.

# Returns
A new `Laser` instance with the given standard deviation, origin, and direction.
"""
Laser(σ, origin::Union{AbstractVector,Tuple}, direction::Union{AbstractVector,Tuple}) = Laser(σ, SVector(Tuple(origin)), SVector(Tuple(direction)))

"""
    sample_source(laser::Laser)

Samples a particle from the given laser source.

# Arguments
- `laser::Laser`: The laser source from which to sample a particle.

# Returns
A particle sampled from the laser source.
"""
function sample_source(laser::Laser)
    radii = √2 * laser.σ * sqrt(randexp())
    θ = 2pi*rand()
    r = @SVector [radii*cos(θ), radii*sin(θ), 0.]

    z_axis = @SVector [0., 0., 1.]
    rot = rotation_a2b(z_axis, laser.direction)

    return Photon(rot*r .+ laser.center, laser.direction)
end