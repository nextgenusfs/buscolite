#!/usr/bin/env python3
"""
BUSCOlite plotting CLI - Generate plots from buscolite JSON results.

This module provides a command-line interface for generating BUSCO assessment
plots from one or more buscolite JSON output files.
"""

import argparse
import json
import os
import sys

from .__init__ import __version__
from .plot import generate_multi_plot, BuscoPlot


def main():
    """
    Main entry point for the buscolite-plot command-line interface.
    
    This function parses command-line arguments and generates plots from
    one or more buscolite JSON result files.
    """
    parser = argparse.ArgumentParser(
        description="Generate BUSCO assessment plots from buscolite JSON results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate plot from a single result
  buscolite-plot sample1.buscolite.json -o sample1_plot.svg
  
  # Generate comparative plot from multiple results
  buscolite-plot sample1.buscolite.json sample2.buscolite.json sample3.buscolite.json -o comparison.svg
  
  # Use wildcards to plot all results in a directory
  buscolite-plot *.buscolite.json -o all_samples.svg
        """,
    )
    
    parser.add_argument(
        "input",
        nargs="+",
        help="One or more buscolite JSON result files",
    )
    
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output SVG file path",
    )
    
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"buscolite-plot v{__version__}",
    )
    
    args = parser.parse_args()
    
    # Validate input files
    for input_file in args.input:
        if not os.path.exists(input_file):
            print(f"Error: Input file not found: {input_file}", file=sys.stderr)
            sys.exit(1)
        if not input_file.endswith(".json"):
            print(f"Warning: Input file does not have .json extension: {input_file}", file=sys.stderr)
    
    # Load JSON data
    datasets = []
    for input_file in args.input:
        try:
            with open(input_file, "r") as f:
                data = json.load(f)
            
            # Extract sample name from filename
            sample_name = os.path.basename(input_file).replace(".buscolite.json", "")
            
            # Check if this is the new format (with stats/config) or old format (just results)
            if "stats" in data and "lineage" in data:
                # New format
                datasets.append({
                    "name": sample_name,
                    "stats": data["stats"],
                    "config": data["lineage"],
                })
            else:
                # Old format - just the results dictionary
                print(f"Error: {input_file} appears to be in old format (missing stats/lineage).", file=sys.stderr)
                print(f"Please re-run buscolite analysis to generate the new JSON format.", file=sys.stderr)
                sys.exit(1)
                
        except json.JSONDecodeError as e:
            print(f"Error: Failed to parse JSON file {input_file}: {e}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Error: Failed to load {input_file}: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Generate plot
    try:
        if len(datasets) == 1:
            # Single dataset - use single plot
            plotter = BuscoPlot(width=800, height=400)
            plotter.create_svg(datasets[0]["stats"], datasets[0]["config"], args.output)
            print(f"✅ Generated single-sample plot: {args.output}")
        else:
            # Multiple datasets - use multi-plot
            generate_multi_plot(datasets, args.output)
            print(f"✅ Generated multi-sample plot: {args.output}")
            print(f"   Compared {len(datasets)} datasets:")
            for ds in datasets:
                stats = ds["stats"]
                complete = stats["single-copy"] + stats["duplicated"]
                total = stats["total"]
                print(f"   - {ds['name']}: C:{complete}/{total} ({complete/total*100:.1f}%)")
    except Exception as e:
        print(f"Error: Failed to generate plot: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

