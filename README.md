# PTMC.jl

[![Build Status](https://github.com/ryarazi/PTMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ryarazi/PTMC.jl/actions/workflows/CI.yml?query=branch%3Amain)

PTMC.jl (Photon Transport Monte Carlo) is a Julia code based on the Monte Carlo modeling of photon transport in multi-layered tissues (MCML) algorithm. The original MCML algorithm was developed by Wang et al. in 1995 with a C code implementation available at [https://omlc.org/software/mc/](https://omlc.org/software/mc/).

The MCML algorithm is a tool for simulating light transport in biological tissues. It is based on the Monte Carlo method, which simulates the path of photons as they travel through the tissue, and is used to calculate the spatial distribution of light fluence, absorption, and scattering in the tissue.

## License

This package is released under the MIT License. Please see the [LICENSE](LICENSE) file for more information.

## Inspiration Sources

A few sources that inspired the development of this package are:

- [MCML: Monte Carlo modeling of light transport in multi-layered tissues](https://omlc.org/software/mc/mcpubs/1995LWCMPBMcml.pdf)
- [Tutorial on Monte Carlo simulation of photon transport in biological tissues](https://opg.optica.org/boe/fulltext.cfm?uri=boe-14-2-559)
- [PyTissueOptics: Python package for Monte Carlo modeling of light transport in biological tissues](https://github.com/DCC-Lab/PyTissueOptics)
