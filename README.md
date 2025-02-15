# misalignedata

This repository contains the code to reproduce the simulation study in the paper "Cohering Disaggregation and Uncertainty Quantification for Spatially Misaligned Data". 

## Abstract


## Code
The code to perform the simulation in the paper consists of two main files, to be executed in the following order in R:

1. `source("compile_mod.R")`: This file loads the necessary libraries, data and compiles the INLA models.
2. `source("score.R")`: This file assesses the models with Squared Error (SE) and Dawid-Sebastiani (DS) scores.

### Description 

1. `load_data.RData`: This file contains the Nepal map and the simulated data.
2. `covariate.R`: This file creates the covariate field.
3. `mesh.R` : This file creates the mesh.
4. `joint_model.R`: This file creates the Observation Plugin (OP) models.
5. `juvpup.R`: This file fits the Joint Uncertainty (JU), Value Plugin (VP) and Uncertainty Plugin (UP) models.
6. `juvpup_nl.R`: This file fits the JU, VP and UP models with non-linear misspecification (NL).
7. `function.R`: This file contains the helper functions used in the paper.