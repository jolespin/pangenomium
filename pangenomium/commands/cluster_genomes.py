"""Cluster genomes into pangenomes using skani"""
from __future__ import print_function, division
import sys
import os
import argparse
from multiprocessing import cpu_count
from collections import OrderedDict
import pandas as pd
from loguru import logger
from pyexeggutor import RunShellCommand, format_header
from .. import __version__
from ..utils import setup_directories, setup_logger, print_header

def parse_input(input_path, genome_extension=None):
    """Parse input as either manifest or simple list, autodetecting format
    
    Accepted organism_types: {prokaryotic, eukaryotic, viral}

    Supports:
    - Manifest (3 cols): [organism_type, id_genome, genome_filepath]
    - Manifest (5+ cols): [organism_type, id_sample, id_genome, genome_filepath, protein_filepath, ...]
    - Simple list: one genome filepath per line
    """
    genome_id_to_filepath = OrderedDict()
    
    # Handle stdin
    if input_path == "stdin" or input_path == "-":
        input_handle = sys.stdin
    else:
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input file not found: {input_path}")
        input_handle = open(input_path, 'r')
    
    try:
        # Read first line to detect format
        first_line = input_handle.readline().strip()
        
        # Check if it's tab-separated (manifest) or not (simple list)
        if '\t' in first_line:
            # It's a manifest - parse with pandas
            if input_path == "stdin" or input_path == "-":
                # Can't seek stdin, so we need to reconstruct
                from io import StringIO
                content = first_line + '\n' + input_handle.read()
                df = pd.read_csv(StringIO(content), sep="\t", header=None)
            else:
                input_handle.close()
                df = pd.read_csv(input_path, sep="\t", header=None)
            
            if df.shape[1] == 3:
                # Batch mode (simple): [organism_type, id_genome, genome_filepath]
                df.columns = ["organism_type", "id_genome", "genome_filepath"]
                for _, row in df.iterrows():
                    genome_id_to_filepath[row["id_genome"]] = row["genome_filepath"]
            
            elif df.shape[1] == 4:
                # Batch mode (full): [organism_type, id_genome, genome_filepath, protein_filepath]
                df.columns = ["organism_type", "id_genome", "genome_filepath", "protein_filepath"]
                for _, row in df.iterrows():
                    genome_id_to_filepath[row["id_genome"]] = row["genome_filepath"]
                    
            elif df.shape[1] >= 5:
                # VEBA mode: [organism_type, id_sample, id_genome, genome_filepath, ...]
                for _, row in df.iterrows():
                    id_genome = row[2]
                    genome_filepath = row[3]
                    genome_id_to_filepath[id_genome] = genome_filepath
            else:
                raise ValueError(f"Manifest must have 3, 4, or ≥5 columns. Got {df.shape[1]}.")
        
        else:
            # It's a simple list - process line by line
            # First pass: collect all filepaths and check extensions
            all_filepaths = [first_line] + [line.strip() for line in input_handle]
            all_filepaths = [fp.strip() for fp in all_filepaths if fp.strip()]
            
            # Validate extension if provided
            if genome_extension:
                ext = genome_extension
                if not ext.startswith("."):
                    ext = "." + ext
                
                # Check how many files have the provided extension
                files_with_ext = [fp for fp in all_filepaths if os.path.basename(fp).endswith(ext)]
                
                if len(files_with_ext) == 0:
                    # No files have this extension - this is an error
                    actual_extensions = set(os.path.splitext(os.path.basename(fp))[1] for fp in all_filepaths)
                    raise ValueError(
                        f"Extension mismatch: You provided extension '{genome_extension}' but none of the "
                        f"{len(all_filepaths)} input files have this extension.\n"
                        f"Actual extensions found: {', '.join(sorted(actual_extensions))}\n"
                        f"Hint: Use '-x {list(actual_extensions)[0][1:]}' for '{list(actual_extensions)[0]}' files"
                    )
                elif len(files_with_ext) < len(all_filepaths):
                    # Some files don't have this extension - warning
                    from loguru import logger
                    logger.warning(
                        f"Extension mismatch: {len(files_with_ext)}/{len(all_filepaths)} files have "
                        f"extension '{ext}'. Files without this extension will be parsed differently."
                    )
            
            # Second pass: process filepaths
            for filepath in all_filepaths:
                basename = os.path.basename(filepath)
                
                if genome_extension:
                    ext = genome_extension
                    if not ext.startswith("."):
                        ext = "." + ext
                    genome_id = basename[:-len(ext)] if basename.endswith(ext) else basename
                else:
                    genome_id = basename.split(".")[0]
                
                genome_id_to_filepath[genome_id] = filepath
    
    finally:
        if input_handle != sys.stdin:
            input_handle.close()
    
    return genome_id_to_filepath

def register_parser(subparsers):
    """Register the cluster-genomes subcommand parser"""
    parser = subparsers.add_parser(
        'cluster-genomes',
        help='Cluster genomes into pangenomes using skani',
        description='Cluster genomes into pangenomes based on Average Nucleotide Identity (ANI)',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # I/O arguments
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument(
        "-i", "--input",
        type=str,
        default="stdin",
        help="Input: manifest file OR list of genome filepaths (one per line) [Default: stdin]"
    )
    parser_io.add_argument("-o", "--output_directory", type=str, default="pangenomium_output/cluster_genomes", help="Output directory [Default: pangenomium_output/cluster_genomes]")
    parser_io.add_argument("-x", "--genome_extension", type=str, help="Genome extension for parsing IDs from list (e.g., 'fa.gz')")
    
    # Utility arguments
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--n_threads", type=int, default=1, help="Number of threads [Default: 1]")
    
    # Skani arguments
    parser_skani = parser.add_argument_group('Skani arguments')
    parser_skani.add_argument("--ani_threshold", type=float, default=95.0, help="ANI threshold [Default: 95.0]")
    parser_skani.add_argument("--minimum_af", type=float, default=50.0, help="Minimum alignment fraction [Default: 50.0]")
    parser_skani.add_argument("--af_mode", type=str, default="relaxed", choices=["relaxed", "strict"], help="AF mode [Default: relaxed]")
    parser_skani.add_argument("--skani_preset", type=str, help="Skani preset")
    parser_skani.add_argument("--skani_options", type=str, default="", help="Additional skani options")
    
    # Clustering arguments
    parser_clustering = parser.add_argument_group('Clustering arguments')
    parser_clustering.add_argument("--organism_type", type=str, choices=["prokaryotic", "eukaryotic", "viral"], help="Organism type (required with --prepend_organism_code)")
    parser_clustering.add_argument("--prepend_organism_code", action="store_true", help="Prepend organism code to cluster prefix (P/E/V for prokaryotic/eukaryotic/viral)")
    parser_clustering.add_argument("--cluster_prefix", type=str, default="SLC-", help="Cluster prefix [Default: 'SLC-']")
    parser_clustering.add_argument("--cluster_suffix", type=str, default="", help="Cluster suffix [Default: '']")
    parser_clustering.add_argument("--cluster_prefix_zfill", type=int, default=0, help="Prefix zfill [Default: 0]")
    parser_clustering.add_argument("--cluster_label_mode", type=str, default="md5", choices=["numeric", "random", "pseudo-random", "md5", "nodes"], help="Label mode [Default: md5]")
    parser_clustering.add_argument("--no_singletons", action="store_true", help="Exclude singletons")
    parser_clustering.add_argument("--identifiers", type=str, help="Identifiers file")
    
    # Set the function to call
    parser.set_defaults(func=run)
    
    return parser

def run(args):
    """Execute cluster-genomes command"""
    
    if args.n_threads == -1:
        args.n_threads = cpu_count()
    
    # Validate organism code usage
    if args.prepend_organism_code and not args.organism_type:
        logger.error("--prepend_organism_code requires --organism_type to be specified")
        return 1
    
    # Prepend organism code to cluster prefix if requested
    cluster_prefix = args.cluster_prefix
    if args.prepend_organism_code and args.organism_type:
        organism_code = args.organism_type[0].upper()  # P, E, or V
        cluster_prefix = organism_code + cluster_prefix
        logger.info(f"Prepending organism code '{organism_code}' to cluster prefix: {cluster_prefix}")
    
    # Setup directories with command-specific subdirectory
    directories = setup_directories(args.output_directory, subdirectory="genome_clustering")
    
    # Create additional output subdirectories
    os.makedirs(os.path.join(directories["output"], "serialization"), exist_ok=True)
    os.makedirs(os.path.join(directories["output"], "representatives"), exist_ok=True)
    
    # Setup logger
    setup_logger(directories["log"], "cluster_genomes.log")
    
    # Print info
    logger.info("="*80)
    logger.info("pangenomium cluster-genomes")
    logger.info("="*80)
    print_header(
        version=__version__,
        n_jobs=args.n_threads,
        additional_info={
            "ANI threshold": args.ani_threshold,
            "Minimum AF": args.minimum_af,
        }
    )
    
    # Parse input
    logger.info("-"*80)
    logger.info("Parsing input")
    logger.info("-"*80)
    input_source = "stdin" if args.input in ["stdin", "-"] else args.input
    logger.info(f"Input source: {input_source}")
    
    genome_id_to_filepath = parse_input(args.input, args.genome_extension)
    logger.info(f"Genomes: {len(genome_id_to_filepath)}")
    
    # Write genome list (filepaths for skani)
    genome_list_filepath = os.path.join(directories["intermediate"], "genome_list.txt")
    with open(genome_list_filepath, "w") as f:
        for genome_id, filepath in genome_id_to_filepath.items():
            print(filepath, file=f)
    
    # Write genome identifiers list (for edgelist-to-clusters)
    genome_identifiers_filepath = os.path.join(directories["intermediate"], "genome_identifiers.list")
    with open(genome_identifiers_filepath, "w") as f:
        for genome_id in genome_id_to_filepath.keys():
            print(genome_id, file=f)
    
    logger.info("")
    
    # ==================
    # Step 1: Run skani
    # ==================
    logger.info("="*80)
    logger.info("Step 1: Running skani triangle")
    logger.info("="*80)
    
    ani_edgelist = os.path.join(directories["intermediate"], "ani_edgelist.tsv")
    
    cmd = [
        "skani", "triangle",
        "-l", genome_list_filepath,
        "-E",
        "-t", str(args.n_threads),
        "-o", ani_edgelist,
    ]
    
    if args.skani_preset:
        cmd.extend(["--preset", args.skani_preset])
    
    if args.skani_options:
        cmd.append(args.skani_options)
    
    # Execute
    step = RunShellCommand(
        command=cmd,
        name="skani",
        validate_output_filepaths=[ani_edgelist]
    ).run()
    step.check_status()
    logger.info("")
    
    # Preprocess skani output: extract first 5 columns and remove header
    logger.info("Preprocessing skani output (extracting columns 1-5)")
    ani_edgelist_processed = os.path.join(directories["intermediate"], "ani_edgelist_processed.tsv")
    
    # Use pandas to extract columns and skip header
    import pandas as pd
    df_ani = pd.read_csv(ani_edgelist, sep="\t", usecols=[0,1,2,3,4])
    df_ani.to_csv(ani_edgelist_processed, sep="\t", index=False, header=False)
    logger.info(f"Processed: {df_ani.shape[0]} edges")
    logger.info("")
    
    # ========================
    # Step 2: Compile clusters
    # ========================
    logger.info("="*80)
    logger.info("Step 2: Compiling genome clusters")
    logger.info("="*80)
    
    genome_clusters = os.path.join(directories["output"], "genomes_to_pangenomes.tsv.gz")
    
    cmd = [
        "edgelist-to-clusters.py",
        "-i", ani_edgelist_processed,
        "-o", genome_clusters,
        "--basename",  # Use basename of paths for matching
        "-t", str(args.ani_threshold),
        "-a", str(args.minimum_af),
        "-m", args.af_mode,
        "--cluster_prefix", cluster_prefix,  # Use cluster_prefix (potentially with organism code prepended)
        "--cluster_prefix_zfill", str(args.cluster_prefix_zfill),
        "--cluster_label_mode", args.cluster_label_mode,
        "--identifiers", genome_identifiers_filepath,  # Provide correct IDs
        "-g", os.path.join(directories["output"], "serialization", "genome_clusters.graph.pkl.gz"),
        "-d", os.path.join(directories["output"], "serialization", "genome_clusters.dict.pkl.gz"),
        "-r", os.path.join(directories["output"], "representatives", "genome_representatives.tsv.gz"),
    ]
    
    # Only add cluster_suffix if non-empty
    if args.cluster_suffix:
        cmd.extend(["--cluster_suffix", args.cluster_suffix])
    
    if args.no_singletons:
        cmd.append("--no_singletons")
    
    # Note: --identifiers is already added above, don't duplicate
    if args.identifiers:
        logger.warning("Ignoring --identifiers argument, using auto-generated identifiers list")
    
    # Execute
    step = RunShellCommand(
        command=cmd,
        name="compile",
        validate_output_filepaths=[genome_clusters]
    ).run()
    step.check_status()
    logger.info("")
    
    logger.info("="*80)
    logger.info("Complete")
    logger.info("="*80)
    return 0
