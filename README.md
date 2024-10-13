[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/buscolite.svg)](https://github.com/nextgenusfs/buscolite/releases/latest)
![Conda](https://img.shields.io/conda/dn/bioconda/buscolite)

# BUSCOlite: simplified BUSCO analysis for gene prediction

BUSCOlite can run the miniprot/Augustus mediated genome predictions as well as the [pyhmmer](https://pyhmmer.readthedocs.io/en/stable/index.html) HMM predictions using the BUSCO v9 or v10 databases. It also provides a python API to run busco analysis from within python, ie to be used inside the eukaryotic gene prediction pipeline Funannotate.

This tool is not meant to be a replacment of BUSCO, for most general use cases you should continue to use [BUSCOv5](https://busco.ezlab.org)

BUSCO models/lineages can be downloaded from the BUSCO site: [v5](https://busco-data.ezlab.org/v5/data/lineages/), [v4](https://busco-data.ezlab.org/v4/data/lineages/).  BUSCOlite does not provide an internal method to do this, as it is trivial to download the lineage you need from your organism(s) by following these links.

##### There are limited dependencies with BUSCOlite:
* [augustus](https://github.com/Gaius-Augustus/Augustus) (note: many versions on conda have non-functional PPX/--proteinprofile mode)
* [miniprot](https://github.com/lh3/miniprot)
* [pyhmmer](https://pyhmmer.readthedocs.io/en/stable/index.html)
* [pyfastx](https://github.com/lmdu/pyfastx)
* [natsort](https://pypi.org/project/natsort/)

#### Why?

[Funannotate](https://github.com/nextgenusfs/funannotate) uses BUSCO to find core conserved marker genes that it uses as a basis to train several ab-initio gene predictors. When BUSCO v2 came out it was python3 only and at that time funannotate was still python2, so I modified the BUSCOv2 source code to be compatible with python2 so it could be run within funannotate. Now BUSCOv5 is the current release, that has numerous bells and whistles that funannotate does not need (no knock against bells and whistles) but the real problem is that due to the large number of dependencies associated with these extra tools is that I cannot build a conda image that includes funannotate and BUSCOv5. So I re-wrote BUSCOv2 here so that it has limited dependencies and will make it easier to incorporate as a dependency of funannotate.  A side note is that the `metaeuk` method that BUSCOv5 now uses as default does not produce complete gene models, in fact the protein sequences it outputs have lowercase sequences that are actually not found in your genome at all.  So for training ab-initio predictors, the `metaeuk` method is not useful -- however, it is faster to get your simple stats on "how complete is my genome assembly".


To install release versions use the pip package manager, like so:
```
python -m pip install buscolite
```

To install the most updated code in master you can run:
```
python -m pip install git+https://github.com/nextgenusfs/buscolite.git
```
