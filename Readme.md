![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Pandas](https://img.shields.io/badge/pandas-%23150458.svg?style=for-the-badge&logo=pandas&logoColor=white)
![NumPy](https://img.shields.io/badge/numpy-%23013243.svg?style=for-the-badge&logo=numpy&logoColor=white)


[![DOI:10.1111/pce.00000](http://img.shields.io/badge/DOI-10.1111/pce.00000-a7d37d.svg)](https://doi.org/10.1111/pce.00000)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# CO Flux article

The code in this repository is complementary to the article "*Unveiling the Influence of Drought and Heat on Leaf Carbon Monoxide Emissions*". The following files are available:

1. **data_full.csv:** Gas flux data of H<sub>2</sub>O, CO<sub>2</sub>, and CO, calculated from concentrations measured using the Aerodyne Mini-TILDAS OCS/COS Monitor QCL Laser ([Link to manual](https://www.aerodyne.com/wp-content/uploads/2021/11/OCS_COS.pdf)). Calculation of fluxes was done according to the [Laser-chamber-fluxes Python scripts](https://github.com/kebasaa/Laser-chamber-fluxes). Contains the following variables:
    - timestamp: Formatted as YYYY-mm-dd HH:MM:SS
	- season: Defined according to [Alpert et al. (2004)](https://doi.org/10.1002/joc.1037)
	- treatment: Droughted & Irrigated
	- co.flux: Carbon monoxide flux, in nmol m<sup>-2</sup> s<sup>-1</sup>
	- co2.flux: CO<sub>2</sub> flux, in umol m<sup>-2</sup> s<sup>-1</sup>
	- Tr: Transpiration flux, in mmol m<sup>-2</sup> s<sup>-1</sup>
	- PAR: Photosynthetically active radiation at the chamber, in umol m<sup>-2</sup> s<sup>-1</sup>
    - PAR_above_canopy: PAR above canopy, in umol m<sup>-2</sup> s<sup>-1</sup>
	- TL: Leaf temperature, in °C
	- TA: Air temperature, in °C
	- VPD: Vapour pressure deficit, in Pandas
	- SWC: Soil water content, averaged from 10-30cm depth, in %
2. **01_CO_figures.ipynb:** Used to create figures for the main article and supplement. Dependencies are:
    - [Pandas](https://pandas.pydata.org/)
    - [Numpy](https://numpy.org/)
    - [Plotnine](https://plotnine.readthedocs.io/en/stable/) & [Mizani](https://plotnine.readthedocs.io/en/stable/tutorials/miscellaneous-manipulating-date-breaks-and-date-labels.html)
3. **02_GAM.R:** Generalised Additive Model development
    - [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html)
	- [ggplot2](https://ggplot2.tidyverse.org/)

## How to Cite

Muller et al. (2023). *Unveiling the Influence of Drought and Heat on Leaf Carbon Monoxide Emissions*. DOI: XYZ

[![DOI:10.1111/pce.00000](http://img.shields.io/badge/DOI-10.1111/pce.00000-a7d37d.svg)](https://doi.org/10.1111/pce.00000)

## License

This software is distributed under the GNU GPL version 3

