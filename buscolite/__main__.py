#!/usr/bin/env python3

import sys
import os
import argparse
from .__version__ import __version__
from .help_formatter import MyParser, MyHelpFormatter
from .busco import runbusco
from .utilities import gffwriter, summary_writer


def main():
    args = parse_args(sys.argv[1:])
    d, m, stats, cfg = runbusco(args.input, args.lineage, species=args.species,
                        mode=args.mode, evalue=args.evalue, cpus=args.cpus,
                        tmpdir='test-buscolite', offset=args.flanks)
    # write gff if genome mode
    if args.mode == 'genome':
        gff = '{}.gff3'.format(args.out)
        with open(gff, 'w') as outfile:
            gffwriter(d, outfile)
    # write summary output
    summary = '{}.buscolite.full_table.tsv'.format(args.out)
    with open(summary, 'w') as outfile:
        summary_writer(d, m, sys.argv, cfg, outfile, mode=args.mode)

    sys.stderr.write('{}\n'.format(stats))




def parse_args(args):
    description = 'BUSCOlite: simplified BUSCO analysis for genome annotation'
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    required.add_argument(
        '-i', '--input', required=True, metavar = '',
        help='Input sequence file in FASTA format (genome or proteome)')
    required.add_argument(
        '-o', '--out', required=True, metavar = '',
        help='Give your analysis run a recognisable short name')
    required.add_argument(
        '-m', '--mode', dest='mode', required=True, metavar = '',
        choices=['genome', 'proteins'],
        help='Specify which BUSCO analysis mode to run. [genome, proteins')
    required.add_argument(
        '-l', '--lineage', required=True, metavar = '',
        help='Specify location of the BUSCO lineage data to be used (full path).')

    optional.add_argument(
        '-c', '--cpus', required=False, type=int, default=1, metavar = '',
        help='Specify the number (N=integer) of threads/cores to use.')
    optional.add_argument(
        '-e', '--evalue', required=False, type=float, default=1e-50, metavar = '',
        help='E-value cutoff for BLAST searches.')
    optional.add_argument(
        '-sp', '--species', required=False, default='anidulans', metavar = '',
        help='Name of existing Augustus species gene finding parameters.')
    optional.add_argument(
        '-f', '--flanks', required=False, type=int, default=2000, metavar = '',
        help='Length of flanking region to use for augustus prediction from tblastn hits.')
    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version',
                            version='{} v{}'.format(
                                os.path.basename(os.path.dirname(os.path.realpath(__file__))),
                                __version__),
                            help="show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


if __name__ == '__main__':
    main()