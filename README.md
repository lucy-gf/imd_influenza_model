# IMD-stratified influenza model

Inferring key influenza-related parameters in England, using a dynamic model 
stratified by Index of Multiple Deprivation (IMD) quintiles, age groups, clinical 
risk groups, and vaccination status. 

---

This project currently uses simulated infection and healthcare data, and tests the 
capacity of the MCMC fitting process to recover the simulation parameters. This 
analysis will then be applied to real-world healthcare data, using the OpenSAFELY 
electronic health records platform.

---

## Overview

The workflow consists of:

1. Producing population data and contact matrices from raw inputs
    - `scripts/setup/load_pop_data.R` and `scripts/setup/load_contact_data.R`
2. Producing "subtype-seasons" data for 2023/24, 2024/25, and 2025/25, using 
[publicly available UKHSA data](https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2025-to-2026/influenza-in-the-uk-annual-epidemiological-report-winter-2025-to-2026#laboratory-surveillance) 
    - `scripts/setup/subtype_setup.R`
3. Producing an `.rds` file of "known" parameters which will be known in the later 
analysis on the OpenSAFELY platform
    - `scripts/produce_known_parameters.R`
    - These include:
        - Vaccination coverage and efficacy against infection and hospitalisation
        - Latent and infectious periods
        - Delays for infections to show up in the healthcare system
        - The proportion of people observed (i.e. in the OpenSAFELY platform)
4. Producing an `.rds` file of "unknown" parameters which will be fitted by MCMC
    - `scripts/produce_unknown_parameters.R`
    - These include:
        - Susceptibility and transmissibility for each subtype-season
        - A parameter controlling epidemic timing (`init_infected`)
        - The proportion of infections each subtype-season in each broad age group, risk group, and 
        IMD quintile who turn up in healthcare systems
5. Simulating three seasons (five subtype-seasons) of influenza infections
    - `scripts/dummy_data/dummy_infections.R`
    - This uses the `run_model()` function in `scripts/seir_model.R`, which is 
    itself using `run_seir_cpp()` (defined in `scripts/seir_model.cpp`)
6. Simulating three seasons of influenza healthcare surveillance data
    - `scripts/dummy_data/dummy_surveillance.R`
    - Healthcare attendances are sampled using a binomial probability (size = 
    number of weekly infections in population subgroup, prob = infection-healthcare-attendance
    ratio of the population subgroup)

> ``📝`` *All of the above can be run using the Makefile, i.e. using the `make` command in the terminal.*

7. Running the MCMC fitting
    - `scripts/dummy_mcmc/mcmc_fitting_HPC.R` or `scripts/dummy_mcmc/mcmc_fitting.R`
        - The output files save differently depending on whether the fitting is 
        run on the MCMC or not; these should otherwise work very similarly
    - This runs the `run_mcmc_inference()` function defined in `scripts/dummy_mcmc/mcmc_functions.R`
    - This is currently set up to take in arguments of `i=1,...,3`, the seasonal index 
    (where `i=1` refers to the 2023/24 season, ..., `i=3` refers to the 2025/26 season)
    and `chain=1,...,10`, running 10 independent chains for each season
       - Each season *jointly* fits to data from the one to two subtypes from that subtype-season

> ``📝`` *The MCMC fitting can be run on the HPC using the bash files `bash/fit_mcmc_i.txt`, where the seasonal index is defined by the `i` in the filename, and the chains are defined as an array.*

8. After running the MCMC fitting on the HPC, extract the outputs (held in `output/data`)
9. Plotting the posterior parameters (traces etc.) 
    - `scripts/dummy_mcmc/plot_mcmc.R`
10. Plotting the posterior epidemics and surveillance data
    - `scripts/dummy_mcmc/plot_mcmc_epids.R`

## Repository Structure

```
imd_influenza_model/
├── data/                    # Input data, with various subfolders
├── scripts/
|   ├── setup/               # Basic functions, color schemes, loading data, etc.
|   ├── dummy_data/          # Produce parameters, infections, surveillance data 
│   ├── dummy_mcmc/          # Run MCMC fitting, analyse posteriors
├── bash/                    # bash files for running MCMC fitting on the HPC
├── mcmc_output/             # .txt files to track MCMC progress
├── renv/                    # R package environment (managed by renv)
├── Makefile                 # Reproducible pipeline (except for fitting contact matrices on HPC)
└── renv.lock                # Locked package versions
```

Output is written to `output/`, with subdirectories for figures (`output/figures/`) and data (`output/data/`).

## Elements to be added

This is work in progress! Examples of elements which are not yet in the model:

``🧓`` Updated age groups (will be changed to 0-4, 5-11, 12-17, 18-29, 30-49, 50-64, 65-74, 75-84, 85+)

``💉`` Better-informed estimates of VE against infection 

``🎆`` Changes in social mixing in the holiday period

``🫄`` Population structure varying in each season

``📊`` Sensitivity analyses around misspecification (e.g. of the vaccine mechanism)

