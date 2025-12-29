"""End-to-end workflow for clustering genomes and proteins"""
from __future__ import print_function, division
import sys
import os
import argparse
from multiprocessing import cpu_count
from loguru import logger
from pyexeggutor import RunShellCommand, format_header
from .. import __version__
from ..utils import setup_directories, setup_logger, print_header

def register_parser(subparsers):
    """Register the end-to-end subcommand parser"""
    parser = subparsers.add_parser(
        'end-to-end',
        help='Run complete pangenomium workflow (genomes + proteins)',
        description='End-to-end workflow: cluster genomes into pangenomes, then cluster proteins within each pangenome',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # I/O arguments
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument(
        "-i", "--input",
        type=str,
        default="stdin",
        help="Input: genome manifest OR list of genome filepaths (one per line) [Default: stdin]"
    )
    parser_io.add_argument("-o", "--output_directory", type=str, default="pangenomium_output", help="Output directory [Default: pangenomium_output]")
    parser_io.add_argument("-x", "--genome_extension", type=str, help="Genome extension for parsing IDs from list (e.g., 'fa.gz')")
    parser_io.add_argument("-m", "--mode", type=str, default="batch", choices=["batch", "veba"], help="Input mode [Default: batch]")
    
    # Utility arguments
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help="Threads [Default: 1]")
    
    # Genome clustering arguments
    parser_genome = parser.add_argument_group('Genome clustering arguments')
    parser_genome.add_argument("--ani_threshold", type=float, default=95.0, help="ANI threshold [Default: 95.0]")
    parser_genome.add_argument("--minimum_af", type=float, default=50.0, help="Minimum AF [Default: 50.0]")
    parser_genome.add_argument("--af_mode", type=str, default="relaxed", choices=["relaxed", "strict"], help="AF mode [Default: relaxed]")
    parser_genome.add_argument("--skani_preset", type=str, help="Skani preset")
    parser_genome.add_argument("--skani_options", type=str, default="", help="Additional skani options")
    parser_genome.add_argument("--genome_cluster_prefix", type=str, default="PanG-", help="Genome cluster prefix [Default: 'PanG-']")
    
    # Protein clustering arguments
    parser_protein = parser.add_argument_group('Protein clustering arguments')
    parser_protein.add_argument("-a", "--algorithm", type=str, default="mmseqs-cluster", choices=["mmseqs-cluster", "mmseqs-linclust"], help="Algorithm [Default: mmseqs-cluster]")
    parser_protein.add_argument("-t", "--minimum_identity_threshold", type=float, default=50.0, help="Identity threshold [Default: 50.0]")
    parser_protein.add_argument("--minimum_coverage_threshold", type=float, default=0.8, help="Coverage threshold [Default: 0.8]")
    parser_protein.add_argument("--mmseqs2_options", type=str, default="", help="MMseqs2 options")
    
    # General clustering arguments
    parser_clustering = parser.add_argument_group('Clustering arguments')
    parser_clustering.add_argument("--cluster_label_mode", type=str, default="md5", choices=["numeric", "random", "pseudo-random", "md5", "nodes"], help="Label mode [Default: md5]")
    parser_clustering.add_argument("--no_singletons", action="store_true", help="Exclude singletons")
    
    # Set the function to call
    parser.set_defaults(func=run)
    
    return parser

def run(args):
    """Execute end-to-end command"""
    
    if args.n_jobs == -1:
        args.n_jobs = cpu_count()
    
    # Setup directories
    directories = setup_directories(args.output_directory)
    
    # Setup logger
    setup_logger(directories["log"], "end_to_end.log")
    
    # Print info
    logger.info("="*80)
    logger.info("pangenomium end-to-end")
    logger.info("="*80)
    print_header(
        version=__version__,
        n_jobs=args.n_jobs,
        additional_info={"Mode": args.mode}
    )
    
    # ===============================
    # Step 1: Run genome clustering
    # ===============================
    logger.info("="*80)
    logger.info("Step 1: Genome clustering")
    logger.info("="*80)
    
    genome_clustering_dir = os.path.join(directories["project"], "genome_clustering")
    
    cmd = [
        "pangenomium", "cluster-genomes",
        "-i", args.input,
        "-o", genome_clustering_dir,
        "-p", str(args.n_jobs),
        "--ani_threshold", str(args.ani_threshold),
        "--minimum_af", str(args.minimum_af),
        "--af_mode", args.af_mode,
        "--cluster_prefix", args.genome_cluster_prefix,
        "--cluster_label_mode", args.cluster_label_mode,
    ]
    
    if args.genome_extension:
        cmd.extend(["-x", args.genome_extension])
    
    if args.skani_preset:
        cmd.extend(["--skani_preset", args.skani_preset])
    
    if args.skani_options:
        cmd.extend(["--skani_options", args.skani_options])
    
    if args.no_singletons:
        cmd.append("--no_singletons")
    
    genome_clusters_file = os.path.join(genome_clustering_dir, "output", "genome_clusters.tsv.gz")
    
    step = RunShellCommand(
        command=cmd,
        name="cluster_genomes",
        validate_output_filepaths=[genome_clusters_file]
    ).run()
    step.check_status()
    logger.info(f"Genome clustering complete: {genome_clusters_file}")
    logger.info("")
    
    # ===============================
    # Step 2: Run protein clustering
    # ===============================
    logger.info("="*80)
    logger.info("Step 2: Protein clustering within pangenomes")
    logger.info("="*80)
    
    protein_clustering_dir = os.path.join(directories["project"], "protein_clustering")
    
    cmd = [
        "pangenomium", "cluster-proteins-from-pangenomes",
        "-g", args.input,
        "-c", genome_clusters_file,
        "-o", protein_clustering_dir,
        "-p", str(args.n_jobs),
        "-t", str(args.minimum_identity_threshold),
        "--minimum_coverage_threshold", str(args.minimum_coverage_threshold),
        "-a", args.algorithm,
        "--cluster_label_mode", args.cluster_label_mode,
    ]
    
    if args.mmseqs2_options:
        cmd.extend(["--mmseqs2_options", args.mmseqs2_options])
    
    if args.no_singletons:
        cmd.append("--no_singletons")
    
    protein_clusters_file = os.path.join(protein_clustering_dir, "output", "protein_clusters.tsv.gz")
    
    step = RunShellCommand(
        command=cmd,
        name="cluster_proteins_from_pangenomes",
        validate_output_filepaths=[protein_clusters_file]
    ).run()
    step.check_status()
    logger.info(f"Protein clustering complete: {protein_clusters_file}")
    logger.info("")
    
    # ==============================
    # Step 3: Create final outputs
    # ==============================
    logger.info("-"*80)
    logger.info("Final outputs")
    logger.info("-"*80)
    
    # Symlink genome clusters
    genome_output = os.path.join(directories["output"], "genome_clusters.tsv.gz")
    if os.path.exists(genome_output):
        os.remove(genome_output)
    os.symlink(os.path.abspath(genome_clusters_file), genome_output)
    logger.info(f"Genome clusters: {genome_output}")
    
    # Symlink protein clusters
    protein_output = os.path.join(directories["output"], "protein_clusters.tsv.gz")
    if os.path.exists(protein_output):
        os.remove(protein_output)
    os.symlink(os.path.abspath(protein_clusters_file), protein_output)
    logger.info(f"Protein clusters: {protein_output}")
    logger.info("")
    
    logger.info("="*80)
    logger.info("Complete")
    logger.info("="*80)
    return 0
