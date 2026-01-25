# BUSCOlite Usage Guide

BUSCOlite provides a command-line interface and a Python API for running BUSCO analysis.

## Command-Line Interface

### Basic Usage

```bash
buscolite -i input.fasta -o output_name -m mode -l lineage_path
```

### Required Arguments

* `-i, --input`: Input sequence file in FASTA format (genome or proteome)
* `-o, --out`: Output name prefix for result files
* `-m, --mode`: Analysis mode, either 'genome' or 'proteins'
* `-l, --lineage`: Path to the BUSCO lineage data

### Optional Arguments

* `-c, --cpus`: Number of CPU threads to use (default: 1)
* `-s, --species`: Augustus species model to use (default: anidulans)
* `-f, --flanks`: Length of flanking regions for gene prediction (default: 2000)
* `-v, --verbose`: Increase verbosity (can be used multiple times)
* `--version`: Show version number and exit
* `-h, --help`: Show help message and exit

### Examples

#### Genome Mode

Analyze a genome using the fungi lineage:

```bash
buscolite -i genome.fasta -o mygenome -m genome -l /path/to/fungi_odb10 -c 8 -s anidulans
```

#### Protein Mode

Analyze a proteome using the fungi lineage:

```bash
buscolite -i proteins.fasta -o myproteins -m proteins -l /path/to/fungi_odb10 -c 8
```

## Output Files

BUSCOlite generates the following output files:

* `<output_name>.buscolite.gff3`: GFF3 file with BUSCO gene annotations (genome mode only)
* `<output_name>.buscolite.tsv`: Tab-separated summary of BUSCO results
* `<output_name>.buscolite.json`: Comprehensive results in JSON format (includes stats, config, and detailed results)

### TSV Output Format

The TSV file contains the following columns:

1. BUSCO ID
2. Status (Complete, Fragmented, Missing)
3. Contig/Sequence ID
4. Start position
5. End position
6. Strand
7. Score
8. Length

### JSON Output Format

The JSON file contains comprehensive information in a structured format:

```json
{
  "version": "25.4.24",
  "mode": "genome",
  "input": "genome.fasta",
  "lineage": {
    "name": "fungi_odb12",
    "creation_date": "2024-01-15",
    "number_of_species": "758",
    "number_of_BUSCOs": "758"
  },
  "stats": {
    "total": 758,
    "single-copy": 720,
    "duplicated": 15,
    "fragmented": 18,
    "missing": 5
  },
  "results": {
    "BUSCO_ID": {
      "status": "complete",
      "contig": "scaffold_1",
      "location": [1000, 2500],
      "strand": "+",
      "hmmer": {
        "bitscore": 450.2,
        "evalue": 1.2e-135
      },
      "translation": "MPROTEINSEQ...",
      ...
    }
  },
  "missing": ["BUSCO_ID1", "BUSCO_ID2"],
  "command": "buscolite -i genome.fasta -l fungi_odb12 -o output"
}
```

This structured format includes:
* **version**: BUSCOlite version used
* **mode**: Analysis mode (genome or proteins)
* **input**: Input file path
* **lineage**: Complete lineage configuration
* **stats**: Summary statistics (total, single-copy, duplicated, fragmented, missing)
* **results**: Detailed information for each BUSCO gene (coordinates, scores, status, sequences, structure)
* **missing**: List of missing BUSCO IDs
* **command**: Command used to generate the results

## Plotting Results

BUSCOlite provides a separate command-line tool for generating publication-quality SVG plots from JSON results.

### buscolite-plot Command

The `buscolite-plot` command can generate plots from one or more BUSCOlite JSON files.

#### Basic Usage

```bash
buscolite-plot <json_file(s)> -o <output.svg>
```

#### Single Sample Plot

Generate a plot from a single analysis:

```bash
buscolite-plot mygenome.buscolite.json -o mygenome_plot.svg
```

This creates a horizontal stacked bar chart showing:
* Complete (Single-copy) - light blue
* Complete (Duplicated) - dark blue
* Fragmented - yellow
* Missing - red

#### Multi-Sample Comparative Plot

Compare multiple samples in a single plot:

```bash
buscolite-plot sample1.buscolite.json sample2.buscolite.json sample3.buscolite.json -o comparison.svg
```

You can also use wildcards to plot all results in a directory:

```bash
buscolite-plot *.buscolite.json -o all_samples.svg
```

The multi-sample plot displays each sample as a separate horizontal bar, making it easy to compare assembly or annotation quality across multiple datasets.

#### Plot Features

* **Zero dependencies**: Uses only Python standard library (SVG generation)
* **Publication quality**: Clean, professional appearance matching BUSCO's official style
* **Automatic detection**: Automatically switches between single and multi-sample modes
* **Informative**: Shows percentages, counts, and summary statistics
* **Lightweight**: Small file sizes, scalable vector graphics

#### Example Workflow

```bash
# Step 1: Run BUSCO analysis on multiple samples
buscolite -i genome1.fasta -o sample1 -m genome -l fungi_odb12 -c 8
buscolite -i genome2.fasta -o sample2 -m genome -l fungi_odb12 -c 8
buscolite -i genome3.fasta -o sample3 -m genome -l fungi_odb12 -c 8

# Step 2: Generate comparative plot
buscolite-plot sample1.buscolite.json sample2.buscolite.json sample3.buscolite.json -o comparison.svg
```

## Python API

BUSCOlite can also be used as a Python library:

```python
from buscolite.busco import runbusco

results, missing, stats, config = runbusco(
    input="genome.fasta",
    lineage="/path/to/fungi_odb10",
    mode="genome",
    species="anidulans",
    cpus=8,
    offset=2000,
    verbosity=3
)

# Print summary
print(f"Complete: {stats['single-copy']}")
print(f"Fragmented: {stats['fragmented']}")
print(f"Missing: {len(missing)}")
print(f"Total: {stats['total']}")

# Access individual BUSCO results
for busco_id, data in results.items():
    if data.get("status") == "complete":
        print(f"{busco_id}: {data['location']}")

# Write results to files
import json
from buscolite.gff import gffwriter
from buscolite.utilities import summary_writer

# Write GFF file
with open("output.gff3", "w") as f:
    gffwriter(results, f)

# Write summary file
with open("output.tsv", "w") as f:
    summary_writer(results, missing, ["test"], config, f, mode="genome")

# Write JSON file
with open("output.json", "w") as f:
    f.write(json.dumps(results, indent=2))
```

For more details on the Python API, see the [API Reference](API.md).
