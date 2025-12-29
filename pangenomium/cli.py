#!/usr/bin/env python
"""Main CLI dispatcher for pangenomium"""
import sys
import argparse
from pangenomium import __version__

def main():
    """Main entry point for pangenomium CLI"""
    
    parser = argparse.ArgumentParser(
        prog='pangenomium',
        description='Dereplicate genomes and proteins into pangenomes and orthologs',
        usage='pangenomium <command> [options]',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'pangenomium v{__version__}'
    )
    
    subparsers = parser.add_subparsers(
        title='commands',
        dest='command',
        help='Available commands',
        metavar='<command>'
    )
    
    # Import command modules
    from pangenomium.commands import (
        cluster_genomes,
        cluster_proteins,
        cluster_proteins_from_pangenomes,
        end_to_end,
        compile_pangenomes_table,
    )
    
    # Register subcommands
    cluster_genomes.register_parser(subparsers)
    cluster_proteins.register_parser(subparsers)
    cluster_proteins_from_pangenomes.register_parser(subparsers)
    end_to_end.register_parser(subparsers)
    compile_pangenomes_table.register_parser(subparsers)
    
    # Parse arguments
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(0)
    
    # Execute the command
    return args.func(args)

if __name__ == '__main__':
    sys.exit(main())
