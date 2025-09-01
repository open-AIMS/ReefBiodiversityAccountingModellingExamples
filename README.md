# Reef Biodiversity Accounting Modelling Examples
Scripts and notebooks demonstrating examples of model runs and analyses for establishing a reef biodiversity accounting project.

## Structure
Each folder coultains a separate `Project.toml` file to activate a Julia environment for running the included scripts.
This can be done by cloning the repo, navigating to the folder you want to run scripts in (e.g. `cd ADRIA_runs`) and activating the project environment via

```
julia --project=.
] instantiate .
```

Each folder also contains Jupyter notebook versions of the scripts with additional explanations and context for the code. For these the environment is activated by running the first cell.

### ADRIA_runs

This folder contains scripts and notebooks  for running intervention scenarios in [ADRIA.jl](https://github.com/open-AIMS/ADRIA.jl) and generating a biodiversity accounting summary table from the results which can be used as input to a Excel/Power BI based analysis.

- `example_ADRIA_runs_biodiversity_accounting_project.jl/ipynb` demonstrates running `ADRIA.jl` for a defined set of intervention scenarios for assessing a biodiversity accoutning project.
- `biodiversity_accounting_table_ADRIA_example.jl/ipynb` demonstrates creating a biodiversity summary table from the set of scenarios run.

### ReefModEngine_runs

This folder contains scripts and notebooks for running intervention scenarios in [ReefModEngine.jl](https://github.com/open-AIMS/ReefModEngine.jl) and generating a biodiversity accounting summary table from the results which can be used as input to a Excel/Power BI based analysis.

- `example_ReefModEngine_runs_biodiversity_accounting_project.jl/ipynb` demonstrates running `ReefModEngine.jl` for a defined set of intervention scenarios for assessing a biodiversity accoutning project.
- `biodiversity_accounting_table_ReefModEngine_example.jl/ipynb` demonstrates creating a biodiversity summary table from the set of scenarios run.
- `plot_ReefModruns_spatial.jl` demonstrates how to cretae spatial plots for a set of `ReefModEngine.jl` results.

### explore_open_access_data

This folder for contains scripts for analysing and visualising open access datasets, including NOAA Reef Watch and Allen Coral Atlas data, to evaluate setting up a reef biodiversity accoutning project. The scripts use [ReefBiodiversityAccountSetup](https://github.com/open-AIMS/ReefBiodiversityAccountSetup.jl) to integrate multiple
spatial data layers for a region of interest, including benthic, geomorphic, depth and historical DHWs.

- `plot_allen_atlas_data_example.ipynb` goes through downloading and plotting spatial datasets from Allen Coral Atlas
and NOAA Reef Watch and combining different data layers into a single dataset of sites identified as suitable for
intervention and control sites. (the `.jl` script does the same, but does not include detailed descriptions of downloading datasets)

- `select_controls_from_ADRIA_data.jl` demonstrates using `ReefBiodiversityAccountSetup` to suggest suitable control sites
for a set of sites selected as intervention sites in a set of `ADRIA.jl` runs.

Note, the functions in `ReefBiodiversityAccountSetup`require a config file detail where Allen Atlas and NOAA data sets are saved. The structure of this `config.toml` file is detailed in `plot_allen_atlas_data_example.ipynb`.

### pluto_notebook_Moore_example
This folder contains a Pluto.jl script for an interactive notebook demonstrting spatial plotting and control and impact site selection for a domain around Moore Reef.

Pluto needs to be installed and then launched in the Julia app to open and run the file. The file also needs the same
`config.toml` file describing data filepaths as decribed in `plot_allen_atlas_data_example.ipynb`.
