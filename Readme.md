[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4117838.svg)](https://doi.org/10.5281/zenodo.4117838)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# CO Flux article

The code in this repository is complementary to the article "Unveiling the Influence of Drought and Heat on Leaf Carbon Monoxide Emissions". The following files are available:

1. data_full.csv: Gas flux data of H<sub>2</sub>O, CO<sub>2</sub>, and CO, calculated from concentrations measured using the Aerodyne Mini-TILDAS OCS/COS Monitor QCL Laser (https://www.aerodyne.com/wp-content/uploads/2021/11/OCS_COS.pdf). Calculation of fluxes was done according to the [Laser-chamber-fluxes Python scripts](https://github.com/kebasaa/Laser-chamber-fluxes).
2. **01_CO_figures.ipynb:** Used to create figures for the main article and supplement. Dependencies are:
  - Pandas
  - Numpy
  - Plotnine & Mizani
3. **02_GAM.R:** Generalised Additive Model development
  - mgcv

## How to Cite

Muller et al. (2023). Unveiling the Influence of Drought and Heat on Leaf Carbon Monoxide Emissions. DOI: 10.5281/zenodo.4117838  (URL:
<https://doi.org/10.5281/zenodo.4117838>), Python notebook

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4117838.svg)](https://doi.org/10.5281/zenodo.4117838)

## License

This software is distributed under the GNU GPL version 3

