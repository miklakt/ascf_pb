# Analytical Self-Consistent Field

A Python package to calculate planar non-charged polymer brush parameters. Such as polymer brush thickness, volume fraction and osmotic pressure at a given distance from grafting surface.

## Abstract

There are next parameters that govern polymer brush's behaviour: grafting density, chain length and affinity to solvent.

At a given distance from grafting surface (z) there is some volume fraction (φ) of the polymer that the brush consists of. One can evaluate this using Analytical Self-Consistent Field approach.

Within the strong stretching approximation there is a universal relation between chemical potential and coordinate z

There is also a normalization condition that says that integrating volume fraction over the brush thickness will give us amount of the polymer in the system.

By using numerical root finding methods it is possible to obtain volume fraction profile.

## Usage

### CLI

The package can be used from CLI with the next syntax:

```
python -m ascf_pb sigma chi N
```

where sigma is grafting density, chi - Flory-Huggins parameter, N - chain length

### Python

After the import you can get a function to evaluate volume fraction or osmotic pressure at a given distance. An instance of such function is initialized with 
sigma, chi and N. For example:

```python
import ascf_pb
phi_func = ascf_pb.phi(sigma=0.05, chi=0.5, N=500)
H = ascf_pb.H(sigma=0.05, chi=0.5, N=500)
z = range(int(H))
phi = phi_func(z)
```

This will calculate phi-profile for the all integer values of distance from zero 
to the brush edge.

## Installation
To install the package run in your python environment
```
pip install [path_to_the_package]
```

## Notes

The package is written as a part of thesis work