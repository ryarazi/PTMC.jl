module PTMC
    using LinearAlgebra
    using Random
    using Rotations
    using StaticArrays
    using MAT
    using Interpolations
    using Plots

    import Plots: plot
    import Base: maximum, minimum, extrema, isapprox, ∈, ∩

    export Laser, Box, Solid, Turbid, Refractive, Scene, HistoryAggregator
    export simulate_particle

    include("linear.jl")

    #geometry
    include("geometry/types.jl")
    include("geometry/ray.jl")
    include("geometry/plane.jl")
    include("geometry/sphere.jl")
    include("geometry/box.jl")
    include("geometry/intersection.jl")

    include("particle.jl")
    include("aggregator.jl")
    include("laser.jl")

    include("material.jl")
    include("solid.jl")
    include("scene.jl")

    include("driver.jl")

    include("plots.jl")
end