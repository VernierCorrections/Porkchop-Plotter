# Porkchop-Plotter
This is a MATLAB script I wrote a year or so ago which can solve Lambert's problem for arbitrary eccentricity values (eg. hyperbolic, parabolic, and elliptical transfers are all supported.
The code itself is not the best, and numerical errors occur in some edge-cases, so I plan on re-writing this in Julia in a more stable manner. 

The included excel sheet includes some orbital parameters for a fictional star system. This is obviously not a very real-world applicable use-case, but demonstrates the code can work for any arrangement (eg. it can be used for transfers between satellites in cislunar orbits, transfers between Jovian moons, interplanetary transfers in our solar system, transfers within the Alpha Centauri system, etc.)

### Callisto-Hops
Callisto-Hops is a modified version of the same code that can use Lambert's problem to solve suborbital transfers, or 'hops', on bodies with little-to-no atmosphere. The code is currently set up with parameters modelling Callisto, one of Jupiter's moons, but could also model our own moon, or Mercury, or any other body with a very thin atmosphere. Suborbital transfers are a very efficient and fast way to get around these bodies, as there is no atmosphere to fly using aerodynamic lift, but conversely there is no drag to worry about.

The attached png files show an example of the output of this script.
