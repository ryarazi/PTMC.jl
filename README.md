# PTMC.jl

[![Build Status](https://github.com/ryarazi/PTMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ryarazi/PTMC.jl/actions/workflows/CI.yml?query=branch%3Amain)

PTMC.jl (Photon Transport Monte Carlo) is a Julia code based on the Monte Carlo modeling of photon transport in multi-layered tissues (MCML) algorithm. The original MCML algorithm was developed by Wang et al. in 1995 with a C code implementation available at [https://omlc.org/software/mc/](https://omlc.org/software/mc/).

The MCML algorithm is a tool for simulating light transport in biological tissues. It is based on the Monte Carlo method, which simulates the path of photons as they travel through the tissue, and is used to calculate the spatial distribution of light fluence, absorption, and scattering in the tissue.

## Example

```julia
using PTMC

# Define the gaussian source parameters
σ = 0.01
direction = [0.0, 0.0, -1.0] #source is pointed in the -̂z direction
center = [0.0, 0.0, 0.02]
source = GuassianSource(σ, direction, center)

# Define geometry
inf_coord = 5000.0 # large enough to be considered infinite
air_box = Box((-inf_coord, -inf_coord, 0.0), (inf_coord, inf_coord, 0.04)) #z from 0 to 0.04
scatterer_box = Box((-inf_coord, -inf_coord, -inf_coord), (inf_coord, inf_coord, 0.0)) #z from -inf (practicaly) to 0

# Define the materials
air_material = Turbid()
scatterer_material = Turbid(μₛ=90, μₐ=10, n=1.5, g=0.0)

# Define the solids
air = Solid(air_box, air_material; name="air")
scatterer = Solid(scatterer_box, scatterer_material; name="scatterer")
scene = Scene([air, scatterer])

#simulate 1000 particles and aggregate their history
aggregator = HistoryAggregator()
for i in 1:1000 simulate_particle(source, scene, aggregator) end

# plot the scene and the history of the first 10 particles
plt = plot(scene)
for i in 1:10
    plot!(aggregator.history[i]; markersize=1.)
end
xlims!(-0.3, 0.3)
ylims!(-0.3, 0.3)
zlims!(-0.05, 0.05)
```

## Inspiration Sources

A few sources that inspired the development of this package are:

- [MCML: Monte Carlo modeling of light transport in multi-layered tissues](https://omlc.org/software/mc/mcpubs/1995LWCMPBMcml.pdf)
- [Tutorial on Monte Carlo simulation of photon transport in biological tissues](https://opg.optica.org/boe/fulltext.cfm?uri=boe-14-2-559)
- [PyTissueOptics: Python package for Monte Carlo modeling of light transport in biological tissues](https://github.com/DCC-Lab/PyTissueOptics)

## License

This package is released under the MIT License. Please see the [LICENSE](LICENSE) file for more information.