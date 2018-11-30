# Microkinetics
A script to solve surface catalytic microkinetics equation 

required models: numpy, scipy, matpoltlib
required files: reactions, species
usage:
python3 ./mk.py

## reactions
reactions file should include all reactions' information.

### notations
* is empty site
A* is adsorbate A on * site
A-B* is transition state
A_g is gas state A

### adsorption and desorption reactions
initial state <-> final state
ex:
* + CO_g <-> CO*
2* + H2_g <-> 2H*
NH3* <-> * + NH3_g

### dissociative adsorption and desorption reactions
initial state <-> transition state <-> final state
ex:
2* + N2_g <-> N-N* + * <-> 2N*
2* + O2_g <-> O-O* + * <-> 2O*

### surface elementary reactions
initial state <-> transition state <-> final state
ex:
N* + H* <-> N-H* + * <-> NH* + *
NH* + H* <-> NH-H* + * <-> NH2* + *
NH2* + H* <-> NH2-H* + * <-> NH3* + *

### gas elementary reactions
initial state <-> transition state <-> final state
ex:
CH4_g + H_g <-> CH4-H_g <-> CH3_g + H2_g
CO_g + O2_g <-> CO-O2_g <-> CO2_g + O_g

## species
species file should include free energies and concentration for all species.
for gas species the concentration is pressure in atm,
for surface species the concentration is initial coverage(Theta).
