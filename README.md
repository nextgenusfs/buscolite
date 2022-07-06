[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/gfftk.svg)](https://github.com/nextgenusfs/buscolite/releases/latest)
![Conda](https://img.shields.io/conda/dn/bioconda/buscolite)

# BUSCOlite: simplified BUSCO analysis for gene prediction

BUSCOlite can run the tblastn/Augustus mediated genome predictions as well as the pyhmmer HMM predictions using the BUSCO v9 or v10 databases. It also provides a python API to run busco analysis from within python, ie to be used inside the eukaryotic gene prediction pipeline Funannotate.

#### Why?

[Funannotate](https://github.com/nextgenusfs/funannotate) uses BUSCO to find core conserved marker genes that it uses as a basis to train several ab-initio gene predictors. When BUSCO v2 came out it was python3 only and at that time funannotate was still python2, so I modified the BUSCOv2 source code to be compatible with python2 so it could be run within funannotate. Now BUSCOv5 is the current release, that has numerous bells and whistles that funannotate does not need (no knock against bells and whistles) but the real problem is that due to the large number of dependencies associated with these extra tools is that I cannot build a conda image that includes funannotate and BUSCOv5. So I re-wrote BUSCOv2 here so that it has limited dependencies and will make it easier to incorporate as a dependency of funannotate.


To install release versions use the pip package manager, like so:
```
python -m pip install buscolite
```

To install the most updated code in master you can run:
```
python -m pip install git+https://github.com/nextgenusfs/buscolite.git
```
