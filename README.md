# Recruitment Analysis
Code and analysis for the paper on juvenile recruitment patterns and their links to adult presence and structural complexity. It stems from VB thesis (Chapter 4).

## Data used
- RR coral annotation from november 2022
- RR DEMs and (orhos) for all sites from november 2022
- settlement tiles from 2022

## Variables and hyphothesis
Adult density: competition or facilitation? (+ or - effects)

Rugosity: hooks and crannies are nice habitats for a coral (+)

Fractal dimension: a fractal surface increases friction and facilitate initial larval settlement (+)

Spatial autocorrelation: patterns are similar for same sites and closer reefs

Settlement tile recruit densities: the more recruit te more juveniles found (+)

## Model and hyphothesis
Gridding reef area into continuous squares of different sizes (50cm, 1m, 2m), computing variables and fitting per each size:

n_juv(nov2022) ~ n_adu(nov2022) + R(nov2022) + D(nov2022) + tile_recruitment(jan2022)+ s(centroid coords) + s(site coords) 

Trying SOAP to accoun for island boundary - SOAP should be well integrated into mgcv (gams), might be trickier to implement in brms but not impossible.

### Ideally also showing 
Realized juvenile settlement areas VS overall means and densities in RR by clipping areas around each juvenile and get R, D, n_adu, n_juv densities and comparing those distributions to overal means.

