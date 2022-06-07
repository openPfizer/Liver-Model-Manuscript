# Liver-Model-Manuscript
QSP model of liver lipid metabolism for 2022 publication.

# Most Recent Release
[![DOI](https://zenodo.org/badge/473639446.svg)](https://zenodo.org/badge/latestdoi/473639446)

# Execution
* `setup.jl` - Will initialize the project, this step can be slow, especially the first time, hence we have split it from the main script.
* `main.jl` - Main project script, executing this script will generate the figures from the manuscript.

# Input Files
* `run_paramaters.csv` - Key parameters for controlling a run of the model, including fit flags and the number of plausible patients to generate
* `parameters_pluto.csv` - Default parameter file, written by executing `derived_parameters.jl` based on initial seeding from `parameters.xlsx`. The values in this file are the ones that used in `main.jl` by default.
* `parameters.xlsx` - Parameter file template that is the input to `derived_parameters.jl` change this file if adding new parameters to the model.

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
