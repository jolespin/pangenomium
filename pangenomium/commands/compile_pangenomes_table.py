"""Compile pangenome manifest table from genome manifest and genome clusters"""
from __future__ import print_function, division
import sys
import os
import argparse
import pandas as pd
from loguru import logger
from .. import __version__

def register_parser(subparsers):
    """Register the compile-pangenomes-table subcommand parser"""
    parser = subparsers.add_parser(
        'compile-pangenomes-table',
        help='Compile pangenome manifest from genome manifest and clusters',
        description="""
Compile pangenome manifest table from genome manifest and genome clusters.

Input formats:
  Batch (4 cols):  [organism_type, id_genome, genome_filepath, protein_filepath]
  Batch (5 cols):  [organism_type, id_genome, genome_filepath, protein_filepath, cds_filepath]
  VEBA (6+ cols):  [organism_type, id_sample, id_genome, genome_filepath, protein_filepath, cds_filepath, ...]
         
Genome clusters: [id_genome, id_pangenome]

Output: pangenome_manifest [id_genome, id_pangenome, protein_filepath]
""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="Genome manifest file (TSV)"
    )
    parser.add_argument(
        "-c", "--clusters",
        type=str,
        required=True,
        help="Genome clusters file (TSV): [id_genome, id_pangenome]"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="stdout",
        help="Output pangenome manifest [Default: stdout]"
    )
    parser.add_argument(
        "--header",
        action="store_true",
        help="Include header in output"
    )
    
    # Set the function to call
    parser.set_defaults(func=run)
    
    return parser

def run(args):
    """Execute compile-pangenomes-table command"""
    
    # Read genome manifest
    df_genomes = pd.read_csv(args.input, sep="\t", header=None)
    
    # Determine format and extract relevant columns
    if df_genomes.shape[1] == 3:
        # Batch format (genomes only): [organism_type, id_genome, genome_filepath]
        # No proteins, so we can't create pangenome manifest
        raise ValueError(
            "Genome manifest has only 3 columns (no protein filepaths). "
            "Pangenome manifest requires protein files. "
            "Use 4-column format: [organism_type, id_genome, genome_filepath, protein_filepath]"
        )
        
    elif df_genomes.shape[1] == 4:
        # Batch format: [organism_type, id_genome, genome_filepath, protein_filepath]
        df_genomes = df_genomes.iloc[:, [1, 3]]
        df_genomes.columns = ["id_genome", "protein_filepath"]
    
    elif df_genomes.shape[1] == 5:
        # Batch format with CDS: [organism_type, id_genome, genome_filepath, protein_filepath, cds_filepath]
        df_genomes = df_genomes.iloc[:, [1, 3]]
        df_genomes.columns = ["id_genome", "protein_filepath"]
        
    elif df_genomes.shape[1] >= 6:
        # VEBA format: [organism_type, id_sample, id_genome, genome_filepath, protein_filepath, ...]
        df_genomes = df_genomes.iloc[:, [2, 4]]
        df_genomes.columns = ["id_genome", "protein_filepath"]
    
    else:
        raise ValueError(
            f"Genome manifest format not recognized. Got {df_genomes.shape[1]} columns. "
            f"Expected: 4-5 columns (batch) or ≥6 columns (VEBA)."
        )
    
    # Read genome clusters
    df_clusters = pd.read_csv(
        args.clusters,
        sep="\t",
        header=None,
        names=["id_genome", "id_pangenome"]
    )
    
    # Merge
    df_output = df_clusters.merge(
        df_genomes,
        on="id_genome",
        how="left"
    )
    
    # Reorder columns: [id_genome, id_pangenome, protein_filepath]
    df_output = df_output[["id_genome", "id_pangenome", "protein_filepath"]]
    
    # Check for missing protein files
    missing = df_output["protein_filepath"].isna().sum()
    if missing > 0:
        logger.warning(f"{missing} genomes have missing protein filepaths")
    
    # Write output
    if args.output == "stdout":
        args.output = sys.stdout
    
    df_output.to_csv(
        args.output,
        sep="\t",
        index=False,
        header=args.header
    )
    
    if args.output != sys.stdout:
        logger.info(f"Pangenome manifest written to: {args.output}")
        logger.info(f"Pangenomes: {df_output['id_pangenome'].nunique()}")
        logger.info(f"Genomes: {len(df_output)}")
    
    return 0
