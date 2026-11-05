
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

# Download the simulation data from [GluA_simulations](https://gin.g-node.org/surbhitw/GluA_simulations). Set `input_folder` in `Scripts/LoadNPlot.py` to the downloaded data location. The dataset identifiers previously used for the three panels are:

```bash
# Figure 5D
python Scripts/LoadNPlot.py -d GluA2/CNIH2_Dep/04_19_2026_18_20_19

# Figure 5E
python Scripts/LoadNPlot.py -d GluA2/CNIH2_Dep/04_19_2026_20_34_10

# Figure 5F
python Scripts/LoadNPlot.py -d  GluA2/CNIH2_Dep/04_19_2026_21_12_22
```

# 4.Contribution

This code was developed by [Surbhit Wagle](https://sites.google.com/view/surbhitwagle/home) with input from all co-authors.

# 4.License

[MIT](https://choosealicense.com/licenses/mit/)

# 5.Data Sources

All data used for plotting are provided in the `Figure_data` folder. Raw data required to rerun some of the analyses can be made available upon reasonable request.

# 6.Acknowledgment

We would like to thank all Tchumathenko, Hafner, Acker-Palmer, Schuman and Choquet lab members for fruitful discussions and helpful feedback.

This study was supported by the University of Bonn Medical Center (SW, CB, NK, TT), the University of Mainz Medical Center (SW, TT), the German Research Foundation via CRC1080 (SW, MKK, AAP, EMS, TT), CRC1507 (AAP), the Donders Institute for Brain, Cognition and Behaviour and Faculty of Science, Radboud University Nijmegen Netherlands (AH). This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (‘MolDynForSyn’, grant agreement No. 945700) (TT) & (‘MemCode’, grant agreement No. 101076961) (AH). AH also received support from the EMBO long-term postdoctoral fellowship (ALTF 1095-2015) and the Alexander von Humboldt Foundation (FRA-1184902-HFST-P).
