This is the read me for the repository containing the scripts used to analyse protein localization and generate figures for the paper
"An integrative model of AMPA receptor trafficking reveals the central contribution of local translation in subtype-specific kinetics"
Surbhit Wagle, Maximilian K. Kracht, Anne Bührke, Amparo Acker-Palmer, Nataliya Kraynyukova, Anne-Sophie Hafner,  Erin M. Schuman,  Tatjana Tchumatchenko.
doi: https://doi.org/10.1101/2025.02.08.637220

# 1. System requirements
The requirement for running this code includes: Python. The code was developed on Apple MacBook Pro M1 running Sonoma 14.5 os with a RAM of 32 GB and
storage of 1 TB. As some of the figure requires downloading simulation data, please insure you have enough RAM and storage space.
# 2.Pre-requisite

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install following packages.

```bash
asteval==0.9.31
Bottleneck @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_29949159-f86f-474b-bc1f-aaa1e0e222b4ofusifik/croots/recipe/bottleneck_1657175564045/work
click @ file:///opt/concourse/worker/volumes/live/2d66025a-4d79-47c4-43be-6220928b6c82/volume/click_1646056610594/work
cloudpickle @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_c57ujq_pgm/croot/cloudpickle_1683040025620/work
contourpy @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_17gskqgptz/croots/recipe/contourpy_1663827415320/work
cycler @ file:///tmp/build/80754af9/cycler_1637851556182/work
cytoolz @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_95v3uo4da9/croot/cytoolz_1667465932724/work
dask @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_ccdxmmxbox/croot/dask-core_1686782920612/work
et-xmlfile==1.1.0
fonttools==4.25.0
fsspec @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_actb8g3z8o/croot/fsspec_1695734513142/work
future==0.18.3
imagecodecs @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_ddpwfftoh7/croot/imagecodecs_1695064957263/work
imageio @ file:///private/var/folders/c_/qfmhj66j0tn016nkx_th4hxm0000gp/T/abs_8994zt3rks/croot/imageio_1695996461372/work
importlib-metadata @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_81_20mq0d8/croot/importlib-metadata_1678997090664/work
kiwisolver @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_e26jwrjf6j/croot/kiwisolver_1672387151391/work
lazy_loader @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_4a85lv11ob/croot/lazy_loader_1695850116141/work
lmfit==1.2.2
locket @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_81b4c56b-0395-495d-81c1-83208d36944d357hqdd0/croots/recipe/locket_1652903116052/work
matplotlib @ file:///private/var/folders/c_/qfmhj66j0tn016nkx_th4hxm0000gp/T/abs_42_ot0zpzy/croot/matplotlib-suite_1693812472014/work
mkl-fft==1.3.6
mkl-random @ file:///Users/ec2-user/mkl/mkl_random_1682994911338/work
mkl-service==2.4.0
munkres==1.1.4
networkx @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_666p3uavvu/croot/networkx_1690562005807/work
numexpr @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_1b50c1js9s/croot/numexpr_1683227065029/work
numpy @ file:///private/var/folders/c_/qfmhj66j0tn016nkx_th4hxm0000gp/T/abs_2ajpp3regc/croot/numpy_and_numpy_base_1691164374110/work
openpyxl==3.0.10
packaging @ file:///private/var/folders/c_/qfmhj66j0tn016nkx_th4hxm0000gp/T/abs_2algm5p9lp/croot/packaging_1693575178038/work
pandas==1.4.2
partd @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_eb38x0gc6u/croot/partd_1693937900739/work
patsy==1.0.1
Pillow==9.4.0
pyparsing @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_3a17y2delq/croots/recipe/pyparsing_1661452538853/work
PyQt5==5.15.9
PyQt5-Qt5==5.15.2
PyQt5-sip==12.12.2
python-dateutil @ file:///tmp/build/80754af9/python-dateutil_1626374649649/work
pytz @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_ddzpsmm2_f/croot/pytz_1671697430473/work
PyWavelets @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_4d73vf63_v/croot/pywavelets_1670425181052/work
PyYAML @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_79xo15pf1i/croot/pyyaml_1670514753622/work
roifile==2023.8.30
scikit-image @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_3fh5dyitqb/croot/scikit-image_1682530834592/work
scikit-posthocs==0.11.1
scipy @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_2eyzr35lpl/croot/scipy_1691606691057/work/dist/scipy-1.11.1-cp310-cp310-macosx_10_9_x86_64.whl#sha256=83f181f25a5e0f2b699e0fb68b754d7aa5bc199091d64fa4c0e3f3621263afd0
seaborn==0.11.2
six @ file:///tmp/build/80754af9/six_1644875935023/work
statannotations==0.5.0
statsmodels==0.14.4
tifffile @ file:///private/var/folders/c_/qfmhj66j0tn016nkx_th4hxm0000gp/T/abs_ffr7rfhtkd/croot/tifffile_1695107463579/work
toolz @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_a7gkswah88/croot/toolz_1667464082910/work
tornado @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_aeotsmw12l/croot/tornado_1690848274212/work
tzdata @ file:///croot/python-tzdata_1690578112552/work
uncertainties==3.1.7
zipp @ file:///private/var/folders/sy/f16zz6x50xz3113nwtb9bvq00000gp/T/abs_b71z79bye2/croot/zipp_1672387125902/work

```
To install the above-listed packages, run: 
```bash
pip install -r requirment.txt 
```
Or you can choose to install them separately

## Usage
Please note that for legacy reasons, we have devided the code base into two repositories.
1. https://github.com/surbhit21/AMPA-dynamics
2. https://github.com/surbhit21/mRNA_analysis.

In order to run the code, you can download or clone this repository, typical installation time is a few minutes (depending on internet speed).
To reproduce the figures, open the repository in terminal and do the following 
(in case you do not see the figure listed below, please check the second repository):

For Figs 2 B and C, S4
```python
python GluA2_protein_analysis.py 
```

For Figs 2 E and F, S5 and S6
```python
python SpineToDend.py
```

For Fig 3 and 5, simulation data needs to be downloaded from https://gin.g-node.org/surbhitw/GluA_simulations.
After that, open LoadNPlot.py file and set up the root_folder (in line 399) to where the downlaoded folder is stored.


For Fig 3 A
```python
python LoadNPlot.py -d 10_25_2024_17_32_15 -s GluA1
```

For Fig 5 D 
```python
python LoadNPlot.py -d 11_08_2024_16_18_55 -s GluA2
```

For Fig 5 E 
```python
python LoadNPlot.py -d 11_11_2024_12_45_55 -s GluA2
```

For Fig 5 F 
```python
python LoadNPlot.py -d 11_11_2024_10_35_01 -s GluA2
```

The run time for each file should take between 10 seconds to a few minutes (depending on your hardware configuration).

# 4.Contribution

This code was developed by [Surbhit Wagle](https://sites.google.com/view/surbhitwagle/home)

# 4.License

[MIT](https://choosealicense.com/licenses/mit/)


# 5.Data Sources
The data to generate figures is provided in the repository. The raw data can be provided upon reasonable request after the peer-review of the manuscript.

# 6.Acknowledgment
This study was supported by the University of Bonn Medical Center (SW, NK, TT), the University of Mainz Medical Center (SW, TT), the German Research Foundation via CRC1080 (SW, TT, MKK, AAP), the Donders Institute for Brain, Cognition and Behaviour and Faculty of Science, Radboud University Nijmegen Netherlands (AH). This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (‘MolDynForSyn’, grant agreement No. 945700) (TT) & (‘MemCode’, grant agreement No. 101076961) (AH). AH also received support from the EMBO long-term postdoctoral fellowship (ALTF 1095-2015) and the Alexander von Humboldt Foundation (FRA-1184902-HFST-P). 


