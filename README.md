
This repository contains the scripts used to run the steady state and time-dependent simulations of GluA1-homomeric and GluA2-containing AMPAR, simulations for the preprint below:

## Preprint

Wagle, S., Kracht, M. K., Bührke, A., Acker-Palmer, A., Kraynyukova, N., Hafner, A.-S., Schuman, E. M., & Tchumatchenko, T. (2025). **An integrative model of AMPA receptor trafficking reveals the central contribution of local translation in subtype-specific kinetics.** *bioRxiv*. https://doi.org/10.1101/2025.02.08.637220

# 1.System requirements

The requirement for running this code includes: Python and jupyter notebook. The code was developed on Apple MacBook Pro M1 running Sonoma 14.5 os with a RAM of 32 GB and
storage of 1 TB. As some of the figure requires downloading simulation data, please insure you have enough RAM and storage space.

# 2.Pre-requisite

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install following packages.

```bash
asteval==0.9.28
contourpy==1.0.6
cycler==0.11.0
fonttools==4.38.0
future==0.18.2
kiwisolver==1.4.4
lmfit==1.0.3
matplotlib==3.6.2
numpy==1.23.4
packaging==21.3
pandas==1.5.1
patsy==0.5.3
Pillow==9.3.0
pyparsing==3.0.9
python-dateutil==2.8.2
pytz==2022.6
scikit-posthocs==0.7.0
scipy==1.9.3
seaborn==0.12.1
six==1.16.0
statsmodels==0.13.5
uncertainties==3.1.7
```

To install the above-listed packages, run:

```bash
pip install -r requirement.txt
```

Or you can choose to install them separately

# 3.Usage

Download or clone this repository and install the dependencies listed above. The table below identifies the analysis entry points for the figure panels. Plotting source data are described in Data Sources; rerunning the original analyses additionally requires the corresponding raw, intermediate, or simulation input files, which are not necessarily included in the repository.

Run Python commands from the repository root (`AMPA-dynamics`), using the `python Scripts/<script>.py` commands below. Python automatically makes `Scripts` available for sibling imports, including `AMPA_model.py`; no `PYTHONPATH` configuration is needed for these commands. Some scripts require input paths and protocol settings to be configured in the source before execution.

| Figure panels | Analysis script or notebook | Analysis | Experimental source data |
| --- | --- | --- | --- |
| 1F–H, S2 | `python Scripts/ParameterExplorationPlots.py` | Effects of transport velocity and surface/internal diffusion coefficients on steady-state AMPAR distributions | Model simulations |
| 3A | `python Scripts/Tanaka2012GluA1plasticity.py`; `compareTanaka()` in `Scripts/LoadNPlot.py` | Simulate GluA1 plasticity and compare with Tanaka and Hirano (2012) | `Source_data/Tanaka_2012_S3A_PSLM.csv` and `Source_data/Tanaka_2012_S3A_NON_PSLM.csv` |
| S4 | `python Scripts/Fitting_to_HaloTagData.py` | HaloTag surface/internal GluA2 profiles, ratios, and model fitting | Profiles in `Data/HaloGluA2_analysis/` |
| 3B, 5D–F, S9, S16 | `python Scripts/GluA2_plasticity.py`; `compareClavet_Fournier()` in `Scripts/LoadNPlot.py` | Generate GluA2 plasticity simulations and compare with Clavet–Fournier PSD data | `Source_data/CF_2023_T4.csv` (Clavet–Fournier 2024 PSD data; see filename note below) |
| S8 | `python Scripts/GluA2_plasticity.py` | Generate GluA2 plasticity simulations using the corresponding protocol settings | Model simulations |
| S7 | `plotTanakaExoprofile()` in `Scripts/GetPublishedData.py` | Plot the Tanaka and Hirano (2012) exocytosis-rate time course and model profile | `Source_data/Tanaka_2012_exo_NPLSM.csv` |

Additional figure-to-script mappings will be added once confirmed.

Use `python Scripts/LoadNPlot.py` to plot saved GluA2 simulation results and compare them with experimental data.

Published experimental data are read by `Scripts/GetPublishedData.py` from `Source_data`, independently of the working directory. Select the appropriate comparison function and simulation dataset in `Scripts/LoadNPlot.py` before running it. For Figure S7, call `plotTanakaExoprofile()`; its example call is currently commented out. Patterson comparisons are not used in the paper.

Filename notes: the Clavet–Fournier PSD table referred to as `CF_2024_T4.csv` is currently stored as `Source_data/CF_2023_T4.csv`. Update the CSV read/write paths in `Scripts/LoadNPlot.py` to that location before plotting; they currently point to `./CF_2023_T4.csv`. The Figure S7 loader expects `Tanaka_2012_exo_NPLSM.csv`, but only `Tanaka_2012_exo_NPLSM_old.csv` is currently present. Confirm the intended data file before generating that panel.

Download the simulation data from [GluA_simulations](https://gin.g-node.org/surbhitw/GluA_simulations). Configure `sim_folder`, `date_time`, `input_folder`, and the protocol metadata path in `Scripts/LoadNPlot.py` for the selected dataset. The current script does not parse a `-d` command-line argument. The dataset identifiers previously used for these panels are:

| Figure panel | Simulation dataset identifier |
| --- | --- |
| 5D | `GluA2/CNIH2_Dep/04_19_2026_18_20_19` |
| 5E | `GluA2/CNIH2_Dep/04_19_2026_20_34_10` |
| 5F | `GluA2/CNIH2_Dep/04_19_2026_21_12_22` |

Use the folder names and capitalization in the downloaded dataset when setting paths.

# 4.Contribution

This code was developed by [Surbhit Wagle](https://sites.google.com/view/surbhitwagle/home) with input from all co-authors.

# 4.License

[MIT](https://choosealicense.com/licenses/mit/)

# 5.Data Sources

Published experimental source data are provided in `Source_data`; processed HaloTag profiles are in `Data/HaloGluA2_analysis`. See the Usage table for figure-specific inputs and filename notes, and the simulation-data link above for saved simulations. Raw data required to rerun some of the analyses can be made available upon reasonable request.

# 6.Acknowledgment

We would like to thank all Tchumathenko, Hafner, Acker-Palmer, Schuman and Choquet lab members for fruitful discussions and helpful feedback.

This study was supported by the University of Bonn Medical Center (SW, CB, NK, TT), the University of Mainz Medical Center (SW, TT), the German Research Foundation via CRC1080 (SW, MKK, AAP, EMS, TT), CRC1507 (AAP), the Donders Institute for Brain, Cognition and Behaviour and Faculty of Science, Radboud University Nijmegen Netherlands (AH). This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (‘MolDynForSyn’, grant agreement No. 945700) (TT) & (‘MemCode’, grant agreement No. 101076961) (AH). AH also received support from the EMBO long-term postdoctoral fellowship (ALTF 1095-2015) and the Alexander von Humboldt Foundation (FRA-1184902-HFST-P).
