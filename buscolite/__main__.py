#!/usr/bin/env python3

import sys
import os
import argparse
import json
from .__init__ import __version__
from .help_formatter import MyParser, MyHelpFormatter
from .busco import runbusco
from .utilities import summary_writer
from .gff import gffwriter
from .log import startLogging


def main():
    args = parse_args(sys.argv[1:])
    logger = startLogging()
    d, m, stats, cfg = runbusco(
        args.input,
        args.lineage,
        species=args.species,
        mode=args.mode,
        cpus=args.cpus,
        offset=args.flanks,
        logger=logger,
        verbosity=3,
    )
    logger.info(
        "Assembly completeness:\n complete={:} [{:.2%}]\n single-copy={:} [{:.2%}]\n fragmented={:} [{:.2%}]\n duplicated={:} [{:.2%}]\n missing={:} [{:.2%}]\n total={:} [{:.2%}]".format(
            stats["single-copy"] + stats["duplicated"],
            ((stats["single-copy"] + stats["duplicated"]) / float(stats["total"])),
            stats["single-copy"],
            (stats["single-copy"] / float(stats["total"])),
            stats["fragmented"],
            (stats["fragmented"] / float(stats["total"])),
            stats["duplicated"],
            (stats["duplicated"] / float(stats["total"])),
            stats["missing"],
            (stats["missing"] / float(stats["total"])),
            stats["total"],
            (stats["total"] / float(stats["total"])),
        )
    )
    # write gff if genome mode
    if args.mode == "genome":
        gff = "{}.buscolite.gff3".format(args.out)
        with open(gff, "w") as outfile:
            gffwriter(d, outfile)
    # write summary output
    summary = "{}.buscolite.tsv".format(args.out)
    with open(summary, "w") as outfile:
        summary_writer(d, m, sys.argv, cfg, outfile, mode=args.mode)
    # write output file to json, might useful
    raw = "{}.buscolite.json".format(args.out)
    with open(raw, "w") as outfile:
        outfile.write(json.dumps(d, indent=2))
    if args.mode == "genome":
        logger.info(
            "Ouput files written:\n GFF3={}\n Summary={}\n Raw={}".format(
                gff, summary, raw
            )
        )
    else:
        logger.info("Ouput files written:\n Summary={}\n Raw={}".format(summary, raw))


def parse_args(args):
    description = "BUSCOlite: simplified BUSCO analysis for genome annotation"
    parser = MyParser(
        description=description, formatter_class=MyHelpFormatter, add_help=False
    )
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")

    required.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="",
        help="Input sequence file in FASTA format (genome or proteome)",
    )
    required.add_argument(
        "-o",
        "--out",
        required=True,
        metavar="",
        help="Give your analysis run a recognisable short name",
    )
    required.add_argument(
        "-m",
        "--mode",
        dest="mode",
        required=True,
        metavar="",
        choices=["genome", "proteins"],
        help="Specify which BUSCO analysis mode to run. [genome, proteins",
    )
    required.add_argument(
        "-l",
        "--lineage",
        required=True,
        metavar="",
        help="Specify location of the BUSCO lineage data to be used (full path).",
    )
    optional.add_argument(
        "-c",
        "--cpus",
        required=False,
        type=int,
        default=1,
        metavar="",
        help="Specify the number (N=integer) of threads/cores to use.",
    )
    optional.add_argument(
        "-s",
        "--species",
        required=False,
        default="anidulans",
        metavar="",
        help="Name of existing Augustus species gene finding parameters.",
    )
    optional.add_argument(
        "-f",
        "--flanks",
        required=False,
        type=int,
        default=2000,
        metavar="",
        help="Length of flanking region to use for augustus prediction from miniprot hits.",
    )
    help_args = parser.add_argument_group("Help")
    help_args.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit",
    )
    help_args.add_argument(
        "--version",
        action="version",
        version="{} v{}".format(
            os.path.basename(os.path.dirname(os.path.realpath(__file__))), __version__
        ),
        help="show program's version number and exit",
    )

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


if __name__ == "__main__":
    main()
