# Analytical Self-Consistent Field

Python package to calculate non-charged polymer brush parameters. 
Such as polymer brush thickness, volume fraction and osmotic pressure 
at given parameters e.g. distance from grafting surface.

## Abstract

There are next parameters that govern polymer brush's behaviour: 
grafting density or number of chains per unit area (sigma), chain length (N) and affinity to solvent (chi). 
Moreover, system geometry(plain, spherical, cylindrical, pore), 
imposed restrictions and brush's topology (linear chains, regular dendron) also influence a polymer brush. 

Within the strong stretching approximation there is a universal relation between 
chemical potential and coordinate z

- log(1 - φ) - 2⋅χ⋅φ = 3/2 κ^2 (Λ^2 - z^2)

There is also a normalization condition that says that integrating volume 
fraction over the brush thickness will give us amount of the polymer 
in the system.

By using numerical root finding methods it is possible to obtain volume 
fraction profile.

## Usage

### CLI

The package can be used from CLI with the next syntax:

```
python -m ascf_pb N sigma chi [-eta eta] [-R R] [-pore_R pore_Radius] [-t topology] [-z_step z_step]
```

positional arguments:
  N                    the chain length of a polymer brush
  sigma                number of chains per unit area
  chi                  Flory-Huggins polymer-solvent interaction parameter

positional arguments:
  N                    the chain length of a polymer brush
  sigma                number of chains per unit area
  chi                  Flory-Huggins polymer-solvent interaction parameter

The default output consist of brush thickness D and
volume fraction profile (an array of z coordinates and an array of volume fraction for a given coordinate)

### Python

After the import you can get a function to evaluate volume fraction or osmotic 
pressure for given parameters. 
The next example shows the case when we want to know how the volume fraction 
changes with varying sigma, having fixed N, chi and z.

```python
import ascf_pb
#generate closure, note we hav not define sigma
phi_from_sigma = ascf_pb.phi(N=1000, chi=0.3, z=10)
#grafting density we want to make calculations with
sigma = [0.01, 0.02, 0.03, 0.04]
#list of resulting local polymer densities
phi_varying_sigma = [phi_from_sigma(s) for s in sigma]
```

This will calculate phi-profile for the all integer values of distance from zero 
to the brush edge.

## Installation
To install the run in your python environment
```
pip install [path_to_the_package]
```

## Notes
The package is written as a part of PhD thesis work of the author.