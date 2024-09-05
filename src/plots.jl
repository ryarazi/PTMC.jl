@recipe function f(b::Box{3,Float64})
    xmin, xmax = minimum(b), maximum(b)
    Δ = xmax - xmin

    xp = [0, 0, 0, 0, 1, 1, 1, 1]*Δ[1] .+ xmin[1]
    yp = [0, 1, 0, 1, 0, 0, 1, 1]*Δ[2] .+ xmin[2]
    zp = [0, 0, 1, 1, 1, 0, 0, 1]*Δ[3] .+ xmin[3]
    connections = [(1,2,3), (4,2,3), (4,7,8), (7,5,6), (2,4,7), (1,6,2), (2,7,6), (7,8,5), (4,8,5), (4,5,3), (1,6,3), (6,3,5)]

    xe = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0]*Δ[1] .+ xmin[1]
    ye = [0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1]*Δ[2] .+ xmin[2]
    ze = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1]*Δ[3] .+ xmin[3]

    xguide --> "x"
    yguide --> "y"
    zguide --> "z"

    @series begin
        label := :none
        linewidth := 0.7
        linecolor := :black
        xe,ye,ze
    end

    @series begin
        fa := 0.3
        linewidth := 0.
        seriestype := :mesh3d
        connections := connections
        xp,yp,zp
    end
end

@recipe function f(s::Solid)
    label --> s.name
    s.geometry
end

@recipe function f(scene::Scene)
    for i in length(scene.solids)-1:-1:1 #first solid is Outside
        @series begin scene.solids[i] end
    end
end

@recipe function f(history::Vector{Photon}; slice = nothing)
    marker --> :circle
    markersize --> 1.
    arrow  --> (:closed, 20.0)
    label --> :none

    @assert slice∈(nothing, "xy", "yx", "yz", "zy", "xz", "zx")

    r = stack(p.r for p in history)
    if isnothing(slice)
        return r[1,:], r[2,:], r[3,:]
    elseif slice=="xy" || slice=="yx"
        return r[1,:], r[2,:]
    elseif slice=="yz" || slice=="zy"
        return r[2,:], r[3,:]
    elseif slice=="zx" || slice=="xz"
        return r[1,:], r[3,:]
    end
end