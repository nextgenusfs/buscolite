# Quick Start

This guide will get you up and running with BUSCOlite in just a few minutes.

## Installation

First, install BUSCOlite:

```bash
pip install buscolite
```

Or using conda:

```bash
conda install -c bioconda buscolite
```

## Download a BUSCO Lineage

Download the appropriate lineage for your organism from the [BUSCO website](https://busco-data.ezlab.org/v5/data/lineages/):

```bash
# Example: Download the fungi lineage
wget https://busco-data.ezlab.org/v5/data/lineages/fungi_odb10.2020-09-10.tar.gz
tar -xzf fungi_odb10.2020-09-10.tar.gz
```

## Run BUSCO Analysis

### Genome Mode

Analyze a genome assembly:

```bash
buscolite -i genome.fasta -o mygenome -m genome -l fungi_odb10 -c 8 -s anidulans
```

This will:

- Analyze `genome.fasta` using the fungi lineage
- Use 8 CPU threads
- Use the `anidulans` Augustus species model
- Generate output files with prefix `mygenome`

### Protein Mode

Analyze a protein set:

```bash
buscolite -i proteins.fasta -o myproteins -m proteins -l fungi_odb10 -c 8
```

## Generate Plots

After running the analysis, you can generate publication-quality SVG plots:

### Single Sample Plot

```bash
buscolite-plot mygenome.buscolite.json -o mygenome_plot.svg
```

### Multi-Sample Comparison

Compare multiple BUSCO results:

```bash
buscolite-plot sample1.buscolite.json sample2.buscolite.json sample3.buscolite.json -o comparison.svg
```

## Understanding the Output

BUSCOlite generates three main output files:

1. **`<output>.buscolite.json`** - Comprehensive results in JSON format including:
   - Version and mode information
   - Lineage configuration
   - Statistics (complete, duplicated, fragmented, missing)
   - Detailed results for each BUSCO
   - Command used

2. **`<output>.buscolite.tsv`** - Tab-separated summary table

3. **`<output>.buscolite.gff3`** - GFF3 annotations (genome mode only)

## Next Steps

- [Usage Guide](usage.md) - Detailed usage documentation
- [API Reference](api/overview.md) - Python API documentation
- [Contributing](contributing.md) - How to contribute to BUSCOlite

