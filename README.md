# Liver-Model-Manuscript
QSP model of liver lipid metabolism for 2022 publication.

# Most Recent Release
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6394196.svg)](https://doi.org/10.5281/zenodo.6394196)

# Execution
* `setup.jl` - Will initialize the project, this step can be slow, especially the first time, hence we have split it from the main script.
* `main.jl` - Main project script, executing this script will generate the figures from the manuscript.

# Additional Functions
* `derived_parameters.jl` - Pluto notebook where we derived some of the baseline parameter values from the model. Load by launching Pluto and opening this notebook.
* `dxdt.jl` - ODEs of the model
* `lit_std_sim.jl` - Script for deriving some values from digitized literature. Not essential for program execution.
* `lognormal_util.jl` - Useful functions for lit_std_sim.jl. Not essential for program execution.
* `mh_sim.jl` - Function for simulating and scoring the model within the Metropolis-Hastings
* `mh.jl` - Main Metropolis-Hastings function used for generating plausible patients.
* `select_vps` - Acceptance/rejection sampling of plausible patients to virtual patients
* `util.jl` - Some shorter, misc. functions that were needed.

![alt text](https://github.com/openPfizer/DigitalHealthData/blob/master/img/osbypfizer.png)
