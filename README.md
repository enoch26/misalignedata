# misalignedata

### Abstract


## Code

### Description 

The code to perform the simulation consists of two main files, to be executed in the following order in R:

1. `source("compile_mod.R")`: This file loads the necessary libraries, data and compiles the INLA models.
2. `source("score.R")`: This file assesses the models with Squared Error (SE) and Dawid-Sebastiani (DS) scores.
