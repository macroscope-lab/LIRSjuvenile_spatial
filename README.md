# Recruitment Analysis
Code and analysis for the paper on juvenile recruitment patterns and their links to adult presence and structural complexity. It stems from VB thesis (Chapter 4).

## Data used
- RR coral annotation for november 2023 and november 2019
- RR DEMs and (orhos) for 15 sites for november 2023 and november 2019
- settlement tiles from jan 2023 and jan 2019

## Variables and hyphothesis
Adult density: competition or facilitation? (+ or - effects) which interacts with other variables

Rugosity: hooks and crannies are nice habitats for a coral (+)

Fractal dimension: a fractal surface increases friction and facilitates initial larval settlement (+)

Settlement tile recruit densities: the more recruits, the more juveniles found (+)

## Model and hypothesis
Gridding reef area into continuous squares of different sizes (50cm, 1m, 2m), computing variables and fitting per each size:

n_juv ~ n_adu*(R+ D+ settlers)+ (1|site) + (1|year)

Tried SOAP to account for island boundary - results hard to interpret and qualitatively not different from site as random effect,

