
# Liao. et al. 2025 Journal of Advances in Modeling Earth Systems

**An unstructured Voronoi mesh framework for large scale hydrologic and land surface models**

Chang Liao<sup>1\*</sup>,
Darren Engwirda<sup>2, 3</sup>,
Donghui Xu<sup>1</sup>,
Tian Zhou<sup>1\*</sup>,
Zeli Tan<sup>1\*</sup>,
and L. Ruby Leung<sup>1</sup>

<sup>1 </sup> Atmospheric Sciences and Global Change, Pacific Northwest National Laboratory, Richland, WA, USA

<sup>2 </sup> T-3 Fluid Dynamics and Solid Mechanics Group, Los Alamos National Laboratory, Los Alamos, NM, USA

<sup>3 </sup> Commonwealth Scientific \& Industrial Research Organisation, Hobart, TAS, Australia

\* corresponding author:  chang.liao@pnnl.gov

## Abstract

Accurate representation of river networks and land surface heterogeneity is crucial for improving the performance of continental-to-global hydrologic and land surface models (LSMs). Traditional approaches rely on either coarse-resolution grids, which may oversimplify complex terrain and limit the accurate representation of key processes, or computationally demanding full-domain high-resolution meshes. This study expands upon our previous research on unstructured mesh-based hydrologic modeling by introducing a novel framework for generating Hydrologic Refined Meshes (HRMs). This framework seamlessly integrates crucial land surface features, including river networks, watershed boundaries, and lake boundaries, into the mesh generation process.
The resulting variable-resolution (4~50km) Voronoi mesh accurately captures these features, demonstrating significant advantages in representing land surface heterogeneity. For instance, the mesh effectively captures surface elevation, resolving most key landscape features despite using nearly half the cell count. River networks are accurately represented by high-resolution cells aligned along their courses, resulting in a conceptual network that deviates by only approximately 1.1km on average from the true river network. We demonstrate the applicability of this framework and evaluate the generated mesh using several key variables, including surface elevation and river networks.

## Journal reference
Liao. et al. (2025). An unstructured Voronoi mesh framework for large scale hydrologic and land surface models.

## Code reference

References for each minted software release for all code involved.

Liao, Chang, & Cooper, Matt. (2022). Pyflowline: a mesh-independent river networks generator for hydrologic models (0.1.22). Zenodo. https://doi.org/10.5281/zenodo.6604337

Liao, Chang. (2022). HexWatershed: a mesh independent flow direction model for hydrologic models (0.1.12). Zenodo. https://doi.org/10.5281/zenodo.6551861


## Data reference

### Input data
Reference for each minted data source for your input data.  For example:



### Output data
Reference for each minted data source for your output data.  For example:



## Contributing modeling software

| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| PyFlowline | version | https://doi.org/10.5281/zenodo.6604337 | 10.5281/zenodo.6604337 |
| HexWatershed | version | https://doi.org/10.5281/zenodo.6551861 | 10.5281/zenodo.6551861 |


## Reproduce my figures

Use the scripts found in the `figures` directory to reproduce the figures used in this publication.


