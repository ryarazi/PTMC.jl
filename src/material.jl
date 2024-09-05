"""
   struct Turbid

A struct to represent a turbid material
"""
struct Turbid
   "material scattering coefficient"
   μₛ::Float64

   "material absorption coefficient"
   μₐ::Float64

   "refractive index"
   n::Float64

   "Henyey-Greenstein coefficient"
   g::Float64

   "total interaction coefficient"
   μₜ::Float64
   
   Turbid(;μₛ=0.,μₐ=0.,n=1.,g=0.) = new(μₛ, μₐ, n, g, μₛ+μₐ)
end

"""
   mean_free_path(mat::Turbid)

Return the mean free path of a turbid material
"""
mean_free_path(mat::Turbid) = mat.μₜ==0. ? Inf : 1/mat.μₜ
Base.show(io::IO, t::Turbid) = print(io, "Turbid(μₛ=$(t.μₛ), μₐ=$(t.μₐ), n=$(t.n), g=$(t.g))")

"""
   Refractive(;n=1.)

Return a refractive material with refractive index `n`
"""
Refractive(;n=1.) = Turbid(n=n)

"""
   Tissue(λₙₘ, B, S, W, M, F, μₛ_reduced_500nm, f, b_mie, g)

Return a turbid material representing a tissue. The model is based on
"Jacques, S. L. (2013). Optical properties of biological tissues: a review". The 
mathematical model is described shortly in  https://omlc.org/news/feb15/generic_optics/index.html.
The code is based on `makeTissueList.m` taken from https://omlc.org/software/mc/mcxyz/

# Arguments
- `λₙₘ::Float64`: wavelength in nm
- `B::Float64`: blood volume fraction
- `S::Float64`: blood oxygen saturation
- `W::Float64`: water volume fraction
- `M::Float64`: melanosome volume fraction
- `F::Float64`: fat volume fraction
- `μₛ_reduced_500nm::Float64`: reduced scattering coefficient at 500nm
- `f::Float64`: scattering power law
- `b_mie::Float64`: scattering power law
- `g::Float64`: Henyey-Greenstein coefficient

# Returns
- `Turbid`: a turbid material representing a tissue
"""

function Tissue(λₙₘ, B, S, W, M, F, μₛ_reduced_500nm, f, b_mie, g)
   data = Dict(label=>vec(mat) for (label,mat) in matread("data/spectralLIB.mat"))
   @assert minimum(data["nmLIB"]) <= λₙₘ <= maximum(data["nmLIB"])

   μₐ_HGb_oxygen = linear_interpolation(data["nmLIB"], data["muaoxy"])(λₙₘ)
   μₐ_HGb_deoxygen = linear_interpolation(data["nmLIB"], data["muadeoxy"])(λₙₘ)
   μₐ_water = linear_interpolation(data["nmLIB"], data["muawater"])(λₙₘ)
   μₐ_melanosome = linear_interpolation(data["nmLIB"], data["muamel"])(λₙₘ)
   μₐ_fat = linear_interpolation(data["nmLIB"], data["muafat"])(λₙₘ)

   μₐ = B*S*μₐ_HGb_oxygen + B*(1-S)*μₐ_HGb_deoxygen + W*μₐ_water + M*μₐ_melanosome + F*μₐ_fat
   μₛ_reduced = μₛ_reduced_500nm * (f*(λₙₘ/500)^(-4) + (1-f)*(λₙₘ/500)^(-b_mie))
   μₛ = μₛ_reduced / (1-g)
   
   ndry = 1.514 #see "Jacques, S. L. (2013). Optical properties of biological tissues: a review"
   nwater = 1.33
   n = ndry - W*(ndry - nwater)

   return Turbid(μₛ=μₛ,μₐ=μₐ,n=n,g=g)
end

#based on makeTissueList.m taken from https://omlc.org/software/mc/mcxyz/
#(see https://omlc.org/software/mc/mcxyz/makeTissueList.m)
Blood(λₙₘ) = Tissue(λₙₘ, 1., 0.75, 0.95, 0., 0., 10., 0., 1.0, 0.9)
Dermis(λₙₘ) = Tissue(λₙₘ, 0.002, 0.67, 0.65, 0., 0., 42.4, 0.62, 1.0, 0.9)
Epidermis(λₙₘ) = Tissue(λₙₘ, 0., 0.75, 0.75, 0.03, 0., 40., 0., 1.0, 0.9)
Skull(λₙₘ) = Tissue(λₙₘ, 0.0005, 0.75, 0.35, 0., 0., 30., 0., 1.0, 0.9)
GrayMatter(λₙₘ) = Tissue(λₙₘ, 0.01, 0.75, 0.75, 0., 0., 20., 0.2, 1.0, 0.9)
WhiteMatter(λₙₘ) = Tissue(λₙₘ, 0.01, 0.75, 0.75, 0., 0., 20., 0.2, 1.0, 0.9)