
"""
    ⟂(x, y)

Check if two vectors are orthogonal.
"""
⟂(x, y) = abs(x⋅y)<1000*eps(eltype(x))

"""
    rotate_vector(u::SVector{3, Float64}, cosθ, ψ)

Rotate a vector `u` by rotaion angles `cosθ` and `ψ`.
"""
function rotate_vector(u::SVector{3, Float64}, cosθ, ψ)
    #@assert norm(u) ≈ 1.
    sinθ = sqrt(1-cosθ^2)
    sinψ, cosψ = sincos(ψ)

    uˣ,uʸ,uᶻ = u
    
    if abs(abs(uᶻ)-1) < 1000*eps(eltype(uᶻ))
        u_rotated = @SVector [sinθ*cosψ, sinθ*sinψ, sign(uᶻ)*cosθ]
    else
        u_rotated = @SVector [sinθ * (uˣ*uᶻ*cosψ - uʸ*sinψ)/sqrt(1-uᶻ^2) + uˣ*cosθ,
                              sinθ * (uʸ*uᶻ*cosψ + uˣ*sinψ)/sqrt(1-uᶻ^2) + uʸ*cosθ,
                             -sinθ * cosψ * sqrt(1-uᶻ^2) + uᶻ*cosθ]
    end
    return u_rotated
end

"""
    rotation_a2b(a::SVector{3, Float64}, b::SVector{3, Float64})

Calculate the rotation matrix that rotates vector `a` to vector `b`.
"""
function rotation_a2b(a::SVector{3, Float64}, b::SVector{3, Float64})
    anorm, bnorm = norm(a), norm(b)
    k = a × b
    cosθ = (a⋅b) / anorm / bnorm
    if abs(cosθ) ≉  1
        v_cross = @SMatrix [0.    -k[3]   k[2]
                            k[3]   0.    -k[1]
                           -k[2]   k[1]   0.]

        return SMatrix{3,3,Float64}(I) + v_cross + (v_cross*v_cross) ./ (1+cosθ) #https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    else # a is parallel to b
        if dot(a,b)>0 return SMatrix{3,3,Float64}(I)
        else return -SMatrix{3,3,Float64}(I) end
    end
end