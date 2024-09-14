"""
    abstract type Source

An abstract type representing a generic source of particles.
"""
abstract type Source end

"""
    struct GuassianSource <: Source

A type representing a gaussian source.

# Fields
- `σ::Float64`: The standard deviation of the source beam.
- `direction::SVector{3,Float64}`: The direction vector of the beam.
- `center::SVector{3,Float64}`: The center position of the gaussian beam.

# Constructor
- `GuassianSource(σ, direction::SVector{3,Float64}, center::SVector{3,Float64})`: Creates a new `GuassianSource` instance with the given standard deviation, direction, and center. The direction vector is normalized.
"""
struct GuassianSource <: Source
    σ::Float64
    direction::SVector{3,Float64}
    center::SVector{3,Float64}

    GuassianSource(σ, direction::SVector{3,Float64}, center::SVector{3,Float64}) = new(σ, normalize(direction), center)
end

"""
    GuassianSource(σ, origin::Union{AbstractVector,Tuple}, direction::Union{AbstractVector,Tuple}) -> GuassianSource

Alternative constructor for `GuassianSource` that accepts the origin and direction as vectors or tuples.

# Arguments
- `σ::Float64`: The standard deviation of the gaussian beam.
- `origin::Union{AbstractVector,Tuple}`: The origin position of the beam.
- `direction::Union{AbstractVector,Tuple}`: The direction vector of the gaussian beam.

# Returns
A new `GuassianSource` instance with the given standard deviation, origin, and direction.
"""
GuassianSource(σ, origin::Union{AbstractVector,Tuple}, direction::Union{AbstractVector,Tuple}) = GuassianSource(σ, SVector(Tuple(origin)), SVector(Tuple(direction)))

"""
    sample_source(source::GuassianSource)

Samples a particle from a given source.

# Arguments
- `source::GuassianSource`: The source from which to sample a particle.

# Returns
A particle sampled from the source.
"""
function sample_source(source::GuassianSource)
    radii = √2 * source.σ * sqrt(randexp())
    θ = 2pi*rand()
    r = @SVector [radii*cos(θ), radii*sin(θ), 0.]

    z_axis = @SVector [0., 0., 1.]
    rot = rotation_a2b(z_axis, source.direction)

    return Photon(rot*r .+ source.center, source.direction)
end