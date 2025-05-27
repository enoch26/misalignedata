# misalignedata

This repository contains the code to reproduce the simulation study in the paper "Coherent Disaggregation and Uncertainty Quantification for Spatially Misaligned Data". 

## Abstract
Spatial misalignment problems arise from both data aggregation and attempts to align misaligned data, leading to information loss. We propose a Bayesian disaggregation framework that links misaligned data to a continuous domain model using an iteratively linearised integration method via integrated nested Laplace approximation (INLA). The framework supports point pattern and aggregated count models under four covariate field scenarios: \textit{Raster at Full Resolution (RastFull), Raster Aggregation (RastAgg), Polygon Aggregation (PolyAgg), and Point Values (PointVal)}. The first three involve aggregation, while the latter two have incomplete fields. For PolyAgg and PointVal, we estimate the full covariate field using \textit{Value Plugin, Joint Uncertainty, and Uncertainty Plugin} methods, with the latter two accounting for uncertainty propagation. These methods demonstrate superior performance, and remain more robust even under model misspecification (i.e. modelling a nonlinear field as linear).

In landslide studies, landslide occurrences are often aggregated into counts based on slope units, reducing spatial detail. The results indicate that point pattern observations and full-resolution covariate fields should be prioritized. For incomplete fields, methods incorporating uncertainty propagation are preferred. This framework supports landslide susceptibility and other spatial mapping, integrating seamlessly with INLA-extension packages.

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

## Citation
For attribution, please cite this work as:
Suen, M. H., Naylor, M., & Lindgren, F. (2025). Coherent Disaggregation and Uncertainty Quantification for Spatially Misaligned Data. arXiv preprint arXiv:2502.10584.

BibTeX citation:
```
@article{suen2025coherent,
  title={Coherent Disaggregation and Uncertainty Quantification for Spatially Misaligned Data},
  author={Suen, Man Ho and Naylor, Mark and Lindgren, Finn},
  journal={arXiv preprint arXiv:2502.10584},
  year={2025}
}
```
