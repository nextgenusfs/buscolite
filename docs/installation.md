# Installation

BUSCOlite can be installed using pip or conda.

## Requirements

BUSCOlite has the following dependencies:

- Python 3.6 or later
- pyhmmer (>=0.12.0)
- pyfastx
- natsort
- augustus (for genome mode)
- miniprot (for genome mode)

## Installing with pip

To install the latest release version using pip:

```bash
pip install buscolite
```

To install the development version directly from GitHub:

```bash
pip install git+https://github.com/nextgenusfs/buscolite.git
```

## Installing with conda

BUSCOlite is available on the bioconda channel:

```bash
conda install -c bioconda buscolite
```

## Installing External Dependencies

For genome mode, you need to install Augustus and miniprot:

### Augustus

Augustus can be installed using conda:

```bash
conda install -c bioconda augustus
```

!!! warning "Augustus PPX Mode"
    Many versions of Augustus on conda have non-functional PPX/--proteinprofile mode. If you encounter issues, you may need to compile Augustus from source.

To compile Augustus from source:

```bash
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus
make
make install
```

### Miniprot

Miniprot can be installed using conda:

```bash
conda install -c bioconda miniprot
```

Or compiled from source:

```bash
git clone https://github.com/lh3/miniprot.git
cd miniprot
make
make install
```

## BUSCO Lineages

BUSCOlite requires BUSCO lineage data, which can be downloaded from the BUSCO website:

- [BUSCO v5 lineages](https://busco-data.ezlab.org/v5/data/lineages/)
- [BUSCO v4 lineages](https://busco-data.ezlab.org/v4/data/lineages/)

Download the appropriate lineage for your organism:

```bash
# Example: Download the fungi lineage
wget https://busco-data.ezlab.org/v5/data/lineages/fungi_odb10.2020-09-10.tar.gz
tar -xzf fungi_odb10.2020-09-10.tar.gz
```

## Verifying Installation

To verify that BUSCOlite is installed correctly:

```bash
buscolite --version
```

This should display the version of BUSCOlite.
