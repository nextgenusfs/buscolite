# API Reference

BUSCOlite provides a Python API for programmatic access to BUSCO analysis functionality.

## Main Modules

### busco

The main module for running BUSCO analysis.

::: buscolite.busco

### search

HMM search functionality using pyhmmer.

::: buscolite.search

### gff

GFF3 file parsing and manipulation.

::: buscolite.gff

### utilities

Utility functions for filtering and processing BUSCO results.

::: buscolite.utilities

### fastx

FASTA/FASTQ file handling.

::: buscolite.fastx

### augustus

Augustus gene prediction integration.

::: buscolite.augustus

### log

Logging utilities.

::: buscolite.log

## Quick Example

```python
from buscolite.busco import runbusco

# Run BUSCO analysis
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
print(f"Complete (single-copy): {stats['single-copy']}")
print(f"Complete (duplicated): {stats['duplicated']}")
print(f"Fragmented: {stats['fragmented']}")
print(f"Missing: {len(missing)}")
print(f"Total: {stats['total']}")

# Access individual BUSCO results
for busco_id, data in results.items():
    if data.get("status") == "complete":
        print(f"{busco_id}: {data['location']}")
```

## Plotting API

```python
from buscolite.plot import generate_plot, generate_multi_plot

# Generate single sample plot
datasets = [{
    'name': 'My Genome',
    'stats': stats,
    'config': config
}]
generate_plot(datasets[0], 'output.svg')

# Generate multi-sample comparison plot
datasets = [
    {'name': 'Sample 1', 'stats': stats1, 'config': config1},
    {'name': 'Sample 2', 'stats': stats2, 'config': config2},
    {'name': 'Sample 3', 'stats': stats3, 'config': config3},
]
generate_multi_plot(datasets, 'comparison.svg')
```

## Detailed Module Documentation

For detailed documentation of each module, see the individual module pages:

- [busco](busco.md) - Main BUSCO analysis functions
- [search](search.md) - HMM search functionality
- [gff](gff.md) - GFF3 file handling
- [utilities](utilities.md) - Utility functions
- [fastx](fastx.md) - FASTA/FASTQ handling
- [augustus](augustus.md) - Augustus integration
- [log](log.md) - Logging utilities
