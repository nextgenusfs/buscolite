# BUSCOlite

**Simplified BUSCO analysis for gene prediction**

BUSCOlite is a streamlined implementation of BUSCO (Benchmarking Universal Single-Copy Orthologs) analysis specifically designed for genome annotation workflows. It can run miniprot/Augustus-mediated genome predictions as well as [pyhmmer](https://pyhmmer.readthedocs.io/en/stable/index.html) HMM predictions using BUSCO v9, v10, or v12 databases.

[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/buscolite.svg)](https://github.com/nextgenusfs/buscolite/releases/latest)
![Conda](https://img.shields.io/conda/dn/bioconda/buscolite)
[![Tests](https://github.com/nextgenusfs/buscolite/actions/workflows/tests.yml/badge.svg)](https://github.com/nextgenusfs/buscolite/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/nextgenusfs/buscolite/branch/master/graph/badge.svg)](https://codecov.io/gh/nextgenusfs/buscolite)

!!! note
    This tool is not meant to be a replacement of BUSCO. For most general use cases you should continue to use [BUSCO](https://busco.ezlab.org).

## Features

- **Genome and protein mode analysis**: Run BUSCO on genome assemblies or protein sets
- **Publication-quality plots**: Generate SVG plots from results with zero additional dependencies
- **Multi-sample comparison**: Compare multiple BUSCO results in a single plot
- **Python API**: Use BUSCOlite programmatically in your own scripts
- **Lightweight**: Minimal dependencies, easy to install and integrate

## Dependencies

BUSCOlite has limited dependencies:

- [augustus](https://github.com/Gaius-Augustus/Augustus) (note: many versions on conda have non-functional PPX/--proteinprofile mode)
- [miniprot](https://github.com/lh3/miniprot)
- [pyhmmer](https://pyhmmer.readthedocs.io/en/stable/index.html)
- [pyfastx](https://github.com/lmdu/pyfastx)
- [natsort](https://pypi.org/project/natsort/)

## Quick Start

Install BUSCOlite:

```bash
pip install buscolite
```

Run BUSCO analysis on a genome:

```bash
buscolite -i genome.fasta -o mygenome -m genome -l /path/to/fungi_odb12 -c 8
```

Generate a plot from the results:

```bash
buscolite-plot mygenome.buscolite.json -o mygenome_plot.svg
```

Compare multiple samples:

```bash
buscolite-plot sample1.buscolite.json sample2.buscolite.json sample3.buscolite.json -o comparison.svg
```

## Why BUSCOlite?

[Funannotate](https://github.com/nextgenusfs/funannotate) uses BUSCO to find core conserved marker genes that it uses as a basis to train several ab-initio gene predictors. When BUSCO v2 came out it was python3 only and at that time funannotate was still python2, so the BUSCO v2 source code was modified to be compatible with python2 so it could be run within funannotate.

Now BUSCOv5 is the current release, which has numerous bells and whistles that funannotate does not need. The real problem is that due to the large number of dependencies associated with these extra tools, it's difficult to build a conda image that includes funannotate and BUSCOv5.

BUSCOlite was re-written to have limited dependencies and make it easier to incorporate as a dependency of funannotate. Additionally, the `metaeuk` method that BUSCOv5 uses as default does not produce complete gene models - the protein sequences it outputs have lowercase sequences that are not actually found in your genome. For training ab-initio predictors, the `metaeuk` method is not useful.

## BUSCO Lineages

BUSCO models/lineages can be downloaded from the BUSCO site:

- [BUSCO v5 lineages](https://busco-data.ezlab.org/v5/data/lineages/)
- [BUSCO v4 lineages](https://busco-data.ezlab.org/v4/data/lineages/)

BUSCOlite does not provide an internal method to download lineages, as it is trivial to download the lineage you need for your organism(s) by following these links.

## Next Steps

- [Installation Guide](installation.md) - Detailed installation instructions
- [Quick Start](quickstart.md) - Get up and running quickly
- [Usage Guide](usage.md) - Comprehensive usage documentation
- [API Reference](api/overview.md) - Python API documentation

## License

BUSCOlite is released under the BSD License. See the [LICENSE](https://github.com/nextgenusfs/buscolite/blob/main/LICENSE.md) file for details.
