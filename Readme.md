![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Pandas](https://img.shields.io/badge/pandas-%23150458.svg?style=for-the-badge&logo=pandas&logoColor=white)
![NumPy](https://img.shields.io/badge/numpy-%23013243.svg?style=for-the-badge&logo=numpy&logoColor=white)


[![DOI:10.1111/nph.20424](http://img.shields.io/badge/DOI-10.1111/nph.20424-a7d37d.svg)](https://doi.org/10.1111/nph.20424)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Data & code for article on leaf CO fluxes

The code in this repository is complementary to the article in New Phytologist titled "*Leaf carbon monoxide emissions under different drought, heat, and light conditions in the field*" (Muller et al., 2025). The following files are available:

1. **data_full.csv:** Gas flux data of H<sub>2</sub>O, CO<sub>2</sub>, and CO, calculated from concentrations measured using the Aerodyne Mini-TILDAS OCS/COS Monitor QCL Laser ([Link to manual](https://www.aerodyne.com/wp-content/uploads/2021/11/OCS_COS.pdf)). Calculation of fluxes was done according to the [Laser-chamber-fluxes Python scripts](https://github.com/kebasaa/Laser-chamber-fluxes). Contains the following variables:
    - _timestamp_: Formatted as YYYY-mm-dd HH:MM:SS
	- _season_: Defined according to [Alpert et al. (2004)](https://doi.org/10.1002/joc.1037)
	- _treatment_: Droughted & Irrigated
	- _co.flux_: Carbon monoxide flux, in nmol m<sup>-2</sup> s<sup>-1</sup>
	- _co2.flux_: CO<sub>2</sub> flux, in μmol m<sup>-2</sup> s<sup>-1</sup>
	- _Tr_: Transpiration flux, in mmol m<sup>-2</sup> s<sup>-1</sup>
	- _PAR_: Photosynthetically active radiation at the chamber, in μmol m<sup>-2</sup> s<sup>-1</sup>
    - _PAR_above_canopy_: PAR above canopy, in μmol m<sup>-2</sup> s<sup>-1</sup>
	- _TL_: Leaf temperature, in °C
	- _TA_: Air temperature, in °C
	- _VPD_: Vapour pressure deficit, in Pa
	- _SWC_: Soil water content, averaged from 10-30cm depth, in %
2. **all_summarised_daily.csv:** Data of *data_full.csv* summarised daily midday
3. **01_CO_figures.ipynb:** Used to create figures for the main article and supplement. Dependencies are:
    - [Pandas](https://pandas.pydata.org/)
    - [Numpy](https://numpy.org/)
    - [Plotnine](https://plotnine.readthedocs.io/en/stable/) & [Mizani](https://plotnine.readthedocs.io/en/stable/tutorials/miscellaneous-manipulating-date-breaks-and-date-labels.html)
4. **02_GAM.R:** Generalised Additive Model development
    - [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html)
	- [ggplot2](https://ggplot2.tidyverse.org/)
	- Other packages, see header of R code file

## How to Cite

Muller et al. (2025). *Leaf carbon monoxide emissions under different drought, heat, and light conditions in the field*. New Phytologist. DOI: https://doi.org/10.1111/nph.20424

[![DOI:10.1111/nph.20424](http://img.shields.io/badge/DOI-10.1111/nph.20424-a7d37d.svg)](https://doi.org/10.1111/nph.20424)

## License

This software is distributed under the GNU GPL version 3

