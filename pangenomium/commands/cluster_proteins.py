"""Cluster proteins into orthogroups using MMseqs2"""
from __future__ import print_function, division
import sys
import os
import argparse
import gzip
from multiprocessing import cpu_count
from loguru import logger
from pyexeggutor import RunShellCommand, format_header
from .. import __version__
from ..utils import setup_directories, setup_logger, print_header

def parse_input(input_path, directories):
    """Parse input as either single FASTA or list of FASTAs, autodetecting format
    
    Supports:
    - Single protein FASTA file
    - List of protein FASTA files (one per line)
    - stdin for piping
    """
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
        
        # Check if first line starts with '>' (FASTA format) or is a filepath
        if first_line.startswith('>'):
            # It's a FASTA file - return the path
            if input_path == "stdin" or input_path == "-":
                # Need to save stdin to a file
                protein_fasta = os.path.join(directories["intermediate"], "stdin_proteins.faa")
                with open(protein_fasta, "w") as f_out:
                    f_out.write(first_line + '\n')
                    f_out.write(input_handle.read())
                return protein_fasta
            else:
                # It's a regular FASTA file
                input_handle.close()
                return input_path
        
        else:
            # It's a list of files - concatenate them
            protein_fasta = os.path.join(directories["intermediate"], "concatenated_proteins.faa")
            
            with open(protein_fasta, "w") as f_out:
                # Process first file
                if first_line:
                    filepath = first_line.strip()
                    if os.path.exists(filepath):
                        if filepath.endswith(".gz"):
                            opener = gzip.open
                            mode = "rt"
                        else:
                            opener = open
                            mode = "r"
                        
                        with opener(filepath, mode) as f_in:
                            for line in f_in:
                                f_out.write(line)
                
                # Process remaining files
                for line in input_handle:
                    filepath = line.strip()
                    if not filepath:
                        continue
                    
                    if not os.path.exists(filepath):
                        logger.warning(f"File not found, skipping: {filepath}")
                        continue
                    
                    if filepath.endswith(".gz"):
                        opener = gzip.open
                        mode = "rt"
                    else:
                        opener = open
                        mode = "r"
                    
                    with opener(filepath, mode) as f_in:
                        for line in f_in:
                            f_out.write(line)
            
            return protein_fasta
    
    finally:
        if input_handle != sys.stdin:
            input_handle.close()

def register_parser(subparsers):
    """Register the cluster-proteins subcommand parser"""
    parser = subparsers.add_parser(
        'cluster-proteins',
        help='Cluster proteins into orthogroups using MMseqs2',
        description='Cluster proteins into orthogroups based on sequence similarity',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # I/O arguments
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument(
        "-i", "--input",
        type=str,
        default="stdin",
        help="Input: protein FASTA file OR list of protein FASTA files (one per line) [Default: stdin]"
    )
    parser_io.add_argument("-o", "--output_directory", type=str, default="pangenomium_output/cluster_proteins", help="Output directory [Default: pangenomium_output/cluster_proteins]")
    
    # Utility arguments
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--n_threads", type=int, default=1, help="Number of threads [Default: 1]")
    
    # MMseqs2 arguments
    parser_mmseqs = parser.add_argument_group('MMseqs2 arguments')
    parser_mmseqs.add_argument("-a", "--algorithm", type=str, default="mmseqs-cluster", choices=["mmseqs-cluster", "mmseqs-linclust"], help="Clustering algorithm [Default: mmseqs-cluster]")
    parser_mmseqs.add_argument("-t", "--minimum_identity_threshold", type=float, default=50.0, help="Identity threshold (0-100) [Default: 50.0]")
    parser_mmseqs.add_argument("-c", "--minimum_coverage_threshold", type=float, default=0.8, help="Coverage threshold (0-1) [Default: 0.8]")
    parser_mmseqs.add_argument("--mmseqs2_options", type=str, default="", help="Additional MMseqs2 options")
    
    # Clustering arguments
    parser_clustering = parser.add_argument_group('Clustering arguments')
    parser_clustering.add_argument("--cluster_prefix", type=str, default="SSPC-", help="Cluster prefix [Default: 'SSPC-']")
    parser_clustering.add_argument("--cluster_suffix", type=str, default="", help="Cluster suffix [Default: '']")
    parser_clustering.add_argument("--cluster_prefix_zfill", type=int, default=0, help="Prefix zfill [Default: 0]")
    parser_clustering.add_argument("--cluster_label_mode", type=str, default="md5", choices=["numeric", "random", "pseudo-random", "md5", "nodes"], help="Label mode [Default: md5]")
    parser_clustering.add_argument("--no_singletons", action="store_true", help="Exclude singletons")
    parser_clustering.add_argument("--identifiers", type=str, help="Identifiers file")
    
    # Set the function to call
    parser.set_defaults(func=run)
    
    return parser

def run(args):
    """Execute cluster-proteins command"""
    
    if args.n_threads == -1:
        args.n_threads = cpu_count()
    
    # Setup directories with command-specific subdirectory
    directories = setup_directories(args.output_directory, subdirectory="protein_clustering")
    
    # Create additional output subdirectories
    os.makedirs(os.path.join(directories["output"], "serialization"), exist_ok=True)
    os.makedirs(os.path.join(directories["output"], "representatives"), exist_ok=True)
    
    # Setup logger
    setup_logger(directories["log"], "cluster_proteins.log")
    
    # Print info
    logger.info("="*80)
    logger.info("pangenomium cluster-proteins")
    logger.info("="*80)
    print_header(
        version=__version__,
        n_jobs=args.n_threads,
        additional_info={
            "Algorithm": args.algorithm,
            "Identity threshold": f"{args.minimum_identity_threshold}%",
            "Coverage threshold": args.minimum_coverage_threshold,
        }
    )
    
    # Parse input
    logger.info("-"*80)
    logger.info("Parsing input")
    logger.info("-"*80)
    input_source = "stdin" if args.input in ["stdin", "-"] else args.input
    logger.info(f"Input source: {input_source}")
    
    protein_fasta = parse_input(args.input, directories)
    logger.info(f"Protein file: {protein_fasta}")
    logger.info("")
    
    # ====================
    # Step 1: Run MMseqs2
    # ====================
    logger.info("="*80)
    logger.info("Step 1: Running MMseqs2")
    logger.info("="*80)
    
    mmseqs_prefix = os.path.join(directories["intermediate"], "mmseqs2")
    protein_edgelist = os.path.join(directories["intermediate"], "protein_edgelist.tsv")
    
    algorithm = "easy-cluster" if args.algorithm == "mmseqs-cluster" else "easy-linclust"
    
    cmd = [
        "mmseqs", algorithm,
        protein_fasta,
        mmseqs_prefix,
        directories["tmp"],
        "--threads", str(args.n_threads),
        "--min-seq-id", str(args.minimum_identity_threshold / 100),
        "-c", str(args.minimum_coverage_threshold),
        "--cov-mode", "1",
    ]
    
    if args.mmseqs2_options:
        cmd.append(args.mmseqs2_options)
    
    cmd.extend([
        "&&",
        "mv", f"{mmseqs_prefix}_cluster.tsv", protein_edgelist,
        "&&",
        "rm", "-rf",
        f"{mmseqs_prefix}_all_seqs.fasta",
        f"{mmseqs_prefix}_rep_seq.fasta",
    ])
    
    # Execute
    step = RunShellCommand(
        command=cmd,
        name="mmseqs2",
        validate_output_filepaths=[protein_edgelist]
    ).run()
    step.check_status()
    logger.info("")
    
    # ========================
    # Step 2: Compile clusters
    # ========================
    logger.info("="*80)
    logger.info("Step 2: Compiling protein clusters")
    logger.info("="*80)
    
    protein_clusters = os.path.join(directories["output"], "protein_clusters.tsv.gz")
    
    cmd = [
        "edgelist-to-clusters.py",
        "-i", protein_edgelist,
        "-o", protein_clusters,
        "--cluster_prefix", args.cluster_prefix,
        "--cluster_prefix_zfill", str(args.cluster_prefix_zfill),
        "--cluster_label_mode", args.cluster_label_mode,
        "-g", os.path.join(directories["output"], "serialization", "protein_clusters.graph.pkl.gz"),
        "-d", os.path.join(directories["output"], "serialization", "protein_clusters.dict.pkl.gz"),
        "-r", os.path.join(directories["output"], "representatives", "protein_representatives.tsv.gz"),
    ]
    
    # Only add cluster_suffix if non-empty
    if args.cluster_suffix:
        cmd.extend(["--cluster_suffix", args.cluster_suffix])
    
    if args.no_singletons:
        cmd.append("--no_singletons")
    
    if args.identifiers:
        cmd.extend(["--identifiers", args.identifiers])
    
    # Execute
    step = RunShellCommand(
        command=cmd,
        name="compile",
        validate_output_filepaths=[protein_clusters]
    ).run()
    step.check_status()
    logger.info("")
    
    logger.info("="*80)
    logger.info("Complete")
    logger.info("="*80)
    return 0
