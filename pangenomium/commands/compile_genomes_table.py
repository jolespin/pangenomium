"""Compile genome manifest table from input directory of genomic assets"""
from __future__ import print_function, division
import sys
import os
import argparse
from collections import defaultdict
import pandas as pd
from loguru import logger
from .. import __version__

def register_parser(subparsers):
    """Register the compile-genomes-table subcommand parser"""
    parser = subparsers.add_parser(
        'compile-genomes-table',
        help='Compile genome manifest from directory of genomic assets',
        description="""
Compile genome manifest table from input directory of genomic assets.

Accepted organism_types: {prokaryotic, eukaryotic, viral}

Output formats:
  Batch mode (4 cols):  [organism_type, id_genome, genome, proteins]
  Batch mode (5 cols):  [organism_type, id_genome, genome, proteins, cds]  (if CDS files detected)
  VEBA mode (7 cols):   [organism_type, id_sample, id_genome, genome, proteins, cds, gff]

Files are discovered by matching basenames: [id_genome].[extension]
Both uncompressed and gzipped files are automatically detected (e.g., .fa and .fa.gz)
""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "-i", "--input_directory",
        type=str,
        required=True,
        help="Input directory containing genomic assets"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="stdout",
        help="Output manifest file [Default: stdout]"
    )
    parser.add_argument(
        "-m", "--mode",
        type=str,
        choices=["batch", "veba"],
        default="batch",
        help="""Input mode:
  batch: [organism_type, id_genome, genome, proteins] (4 cols) OR
         [organism_type, id_genome, genome, proteins, cds] (5 cols)
  veba:  [organism_type, id_sample, id_genome, genome, proteins, cds, gff] (7 cols)
[Default: batch]"""
    )
    
    # Mutually exclusive organism type specification
    organism_group = parser.add_mutually_exclusive_group(required=True)
    organism_group.add_argument(
        "--genome_to_organism",
        type=str,
        help="TSV file mapping [id_genome, organism_type] (no header)"
    )
    organism_group.add_argument(
        "-t", "--organism_type",
        type=str,
        choices=["prokaryotic", "eukaryotic", "viral"],
        help="Organism type for all genomes (mutually exclusive with --genome_to_organism)"
    )
    
    # VEBA mode requirement
    parser.add_argument(
        "--genome_to_sample",
        type=str,
        help="TSV file mapping [id_genome, id_sample] (no header) - Required for VEBA mode"
    )
    
    # File extensions
    parser.add_argument(
        "-a", "--assembly_extension",
        type=str,
        default=".fa",
        help="Assembly/genome file extension (also detects .gz) [Default: .fa]"
    )
    parser.add_argument(
        "-p", "--protein_extension",
        type=str,
        default=".faa",
        help="Protein file extension (also detects .gz) [Default: .faa]"
    )
    parser.add_argument(
        "-c", "--cds_extension",
        type=str,
        default=".ffn",
        help="CDS file extension (also detects .gz) [Default: .ffn]"
    )
    parser.add_argument(
        "-g", "--gff_extension",
        type=str,
        default=".gff",
        help="GFF file extension (also detects .gz) [Default: .gff]"
    )
    
    # Set the function to call
    parser.set_defaults(func=run)
    
    return parser

def run(opts):
    """Execute compile-genomes-table command"""
    
    # Validate mode-specific requirements
    if opts.mode == "veba" and not opts.genome_to_sample:
        logger.error("VEBA mode requires --genome_to_sample mapping file")
        return 1
    
    # Validate input directory
    if not os.path.isdir(opts.input_directory):
        logger.error(f"Input directory does not exist: {opts.input_directory}")
        return 1
    
    logger.info(f"Scanning directory: {opts.input_directory}")
    logger.info(f"Mode: {opts.mode}")
    
    # Helper function to check both uncompressed and gzipped versions
    def check_extension(fname, ext):
        """Check if filename matches extension or extension.gz and return genome_id"""
        if fname.endswith(ext + ".gz"):
            return fname[:-len(ext + ".gz")]
        elif fname.endswith(ext):
            return fname[:-len(ext)]
        return None
    
    # Scan directory and group files by genome ID
    genome_files = defaultdict(dict)
    
    for filename in os.listdir(opts.input_directory):
        filepath = os.path.join(opts.input_directory, filename)
        
        # Skip directories
        if os.path.isdir(filepath):
            continue
        
        # Extract genome ID by removing extension (handles .gz automatically)
        genome_id = check_extension(filename, opts.assembly_extension)
        if genome_id:
            genome_files[genome_id]['genome'] = filepath
            continue
            
        genome_id = check_extension(filename, opts.protein_extension)
        if genome_id:
            genome_files[genome_id]['proteins'] = filepath
            continue
            
        genome_id = check_extension(filename, opts.cds_extension)
        if genome_id:
            genome_files[genome_id]['cds'] = filepath
            continue
            
        genome_id = check_extension(filename, opts.gff_extension)
        if genome_id:
            genome_files[genome_id]['gff'] = filepath
    
    logger.info(f"Discovered {len(genome_files)} unique genome IDs")
    
    # Debug: Log detection summary
    if genome_files:
        sample_genome = list(genome_files.keys())[0]
        logger.info(f"Example genome '{sample_genome}' has files: {list(genome_files[sample_genome].keys())}")
        
        file_type_counts = {
            'genome': sum(1 for f in genome_files.values() if 'genome' in f),
            'proteins': sum(1 for f in genome_files.values() if 'proteins' in f),
            'cds': sum(1 for f in genome_files.values() if 'cds' in f),
            'gff': sum(1 for f in genome_files.values() if 'gff' in f),
        }
        logger.info(f"File type counts: genome={file_type_counts['genome']}, proteins={file_type_counts['proteins']}, cds={file_type_counts['cds']}, gff={file_type_counts['gff']}")
    
    # Load organism type mapping
    if opts.genome_to_organism:
        logger.info(f"Loading organism type mapping: {opts.genome_to_organism}")
        df_organism = pd.read_csv(
            opts.genome_to_organism,
            sep="\t",
            header=None,
            names=["id_genome", "organism_type"]
        )
        organism_map = dict(zip(df_organism['id_genome'], df_organism['organism_type']))
    else:
        logger.info(f"Using organism type: {opts.organism_type}")
        organism_map = {genome_id: opts.organism_type for genome_id in genome_files.keys()}
    
    # Load sample mapping for VEBA mode
    if opts.mode == "veba":
        logger.info(f"Loading genome-to-sample mapping: {opts.genome_to_sample}")
        df_sample = pd.read_csv(
            opts.genome_to_sample,
            sep="\t",
            header=None,
            names=["id_genome", "id_sample"]
        )
        sample_map = dict(zip(df_sample['id_genome'], df_sample['id_sample']))
        
        # Validate all genomes have sample mapping
        missing_samples = set(genome_files.keys()) - set(sample_map.keys())
        if missing_samples:
            logger.error(f"{len(missing_samples)} genomes missing sample mapping")
            logger.error(f"Example missing genomes: {list(missing_samples)[:5]}")
            return 1
    
    # Build manifest data
    rows = []
    missing_genomes = []
    missing_proteins = []
    missing_cds = []
    missing_gff = []
    
    # Check if any genome has CDS files (for batch mode column decision)
    has_cds = any('cds' in files for files in genome_files.values())
    
    if opts.mode == "batch":
        logger.info(f"CDS files detected: {has_cds}")
        # In batch mode, if ANY genome has CDS, ALL must have CDS (otherwise inconsistent output)
        if has_cds:
            genomes_without_cds = [gid for gid, files in genome_files.items() if 'cds' not in files and 'genome' in files and 'proteins' in files]
            if genomes_without_cds:
                logger.error(f"Batch mode detected CDS files, but {len(genomes_without_cds)} genomes are missing CDS files")
                logger.error(f"Either all genomes must have CDS or none should have CDS")
                logger.error(f"Example genomes without CDS: {genomes_without_cds[:5]}")
                return 1
    
    for genome_id, files in sorted(genome_files.items()):
        # Validate required files
        if 'genome' not in files:
            missing_genomes.append(genome_id)
            continue
        if 'proteins' not in files:
            missing_proteins.append(genome_id)
            continue
        
        # VEBA mode requires CDS and GFF
        if opts.mode == "veba":
            if 'cds' not in files:
                missing_cds.append(genome_id)
                continue
            if 'gff' not in files:
                missing_gff.append(genome_id)
                continue
        
        # Get organism type
        organism_type = organism_map.get(genome_id, "unknown")
        
        if opts.mode == "batch":
            if has_cds:
                # 5-column format
                row = [
                    organism_type,
                    genome_id,
                    files['genome'],
                    files['proteins'],
                    files['cds']  # CDS is required if has_cds is True
                ]
            else:
                # 4-column format
                row = [
                    organism_type,
                    genome_id,
                    files['genome'],
                    files['proteins']
                ]
        else:  # VEBA mode
            # 7-column format - all files are required in VEBA mode
            row = [
                organism_type,
                sample_map[genome_id],
                genome_id,
                files['genome'],
                files['proteins'],
                files['cds'],
                files['gff']
            ]
        
        rows.append(row)
    
    # Report missing files
    if missing_genomes:
        logger.error(f"{len(missing_genomes)} genomes missing assembly files")
        return 1
    if missing_proteins:
        logger.error(f"{len(missing_proteins)} genomes missing protein files")
        return 1
    
    if opts.mode == "veba":
        if missing_cds:
            logger.error(f"{len(missing_cds)} genomes missing CDS files")
            return 1
        if missing_gff:
            logger.error(f"{len(missing_gff)} genomes missing GFF files")
            return 1
    
    # Create DataFrame
    df_output = pd.DataFrame(rows)
    
    # Write output without header or index
    if opts.output == "stdout":
        output_handle = sys.stdout
    else:
        output_handle = opts.output
    
    df_output.to_csv(
        output_handle,
        sep="\t",
        index=False,
        header=False
    )
    
    if opts.output != "stdout":
        logger.info(f"Genome manifest written to: {opts.output}")
        logger.info(f"Total genomes: {len(df_output)}")
        if opts.mode == "batch":
            logger.info(f"Format: {len(df_output.columns)} columns")
        else:
            logger.info(f"Samples: {df_output.iloc[:, 1].nunique()}")
    
    return 0
