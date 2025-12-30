"""Cluster proteins within pangenomes into pangenome-specific orthogroups"""
from __future__ import print_function, division
import sys
import os
import argparse
import gzip
from multiprocessing import cpu_count
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
from loguru import logger
from pyexeggutor import RunShellCommand, format_header
from .. import __version__
from ..utils import setup_directories, setup_logger, print_header

def parse_input(args, directories):
    """Parse pangenome manifest and prepare protein files"""
    pangenome_to_proteins = defaultdict(dict)
    
    if args.pangenome_manifest:
        # Direct pangenome manifest: [id_genome, id_pangenome, protein_filepath]
        df = pd.read_csv(args.pangenome_manifest, sep="\t", header=None)
        df.columns = ["id_genome", "id_pangenome", "protein_filepath"]
        
        for pangenome_id, group in df.groupby("id_pangenome"):
            # Create concatenated protein file for this pangenome
            protein_fasta = os.path.join(
                directories["intermediate"],
                f"{pangenome_id}_proteins.faa"
            )
            
            with open(protein_fasta, "w") as f_out:
                for _, row in group.iterrows():
                    filepath = row["protein_filepath"]
                    
                    if filepath.endswith(".gz"):
                        opener = gzip.open
                        mode = "rt"
                    else:
                        opener = open
                        mode = "r"
                    
                    with opener(filepath, mode) as f_in:
                        for line in f_in:
                            f_out.write(line)
            
            pangenome_to_proteins[pangenome_id]["fasta"] = protein_fasta
            pangenome_to_proteins[pangenome_id]["genomes"] = group["id_genome"].tolist()
    
    elif args.genome_manifest and args.genome_clusters:
        # Create pangenome manifest from genome manifest + clusters
        df_genomes = pd.read_csv(args.genome_manifest, sep="\t", header=None)
        df_clusters = pd.read_csv(args.genome_clusters, sep="\t", header=None, names=["id_genome", "id_pangenome"])
        
        # Determine manifest format
        if df_genomes.shape[1] == 4:
            # Batch: [organism_type, id_genome, genome_filepath, protein_filepath]
            df_genomes.columns = ["organism_type", "id_genome", "genome_filepath", "protein_filepath"]
            df_genomes = df_genomes[["id_genome", "protein_filepath"]]
        elif df_genomes.shape[1] >= 5:
            # VEBA: [organism_type, id_sample, id_genome, genome_filepath, protein_filepath, ...]
            df_genomes = df_genomes.iloc[:, [2, 4]]
            df_genomes.columns = ["id_genome", "protein_filepath"]
        
        # Merge with clusters
        df_merged = df_clusters.merge(df_genomes, on="id_genome", how="left")
        
        for pangenome_id, group in df_merged.groupby("id_pangenome"):
            protein_fasta = os.path.join(
                directories["intermediate"],
                f"{pangenome_id}_proteins.faa"
            )
            
            with open(protein_fasta, "w") as f_out:
                for _, row in group.iterrows():
                    filepath = row["protein_filepath"]
                    
                    if pd.isna(filepath) or not os.path.exists(filepath):
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
            
            pangenome_to_proteins[pangenome_id]["fasta"] = protein_fasta
            pangenome_to_proteins[pangenome_id]["genomes"] = group["id_genome"].tolist()
    
    return pangenome_to_proteins

def register_parser(subparsers):
    """Register the cluster-proteins-from-pangenomes subcommand parser"""
    parser = subparsers.add_parser(
        'cluster-proteins-from-pangenomes',
        help='Cluster proteins within pangenomes into pangenome-specific orthogroups',
        description='Cluster proteins within each pangenome separately using MMseqs2',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # I/O arguments
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-i", "--pangenome_manifest", type=str, help="Pangenome manifest: [id_genome, id_pangenome, protein_filepath]")
    parser_io.add_argument("-g", "--genome_manifest", type=str, help="Genome manifest (use with --genome_clusters)")
    parser_io.add_argument("-c", "--genome_clusters", type=str, help="Genome clusters file (use with --genome_manifest)")
    parser_io.add_argument("-o", "--output_directory", type=str, default="pangenome_protein_clustering_output", help="Output directory [Default: pangenome_protein_clustering_output]")
    
    # Utility arguments
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("-p", "--n_jobs", type=int, default=1, help="Threads [Default: 1]")
    
    # MMseqs2 arguments
    parser_mmseqs = parser.add_argument_group('MMseqs2 arguments')
    parser_mmseqs.add_argument("-a", "--algorithm", type=str, default="mmseqs-cluster", choices=["mmseqs-cluster", "mmseqs-linclust"], help="Algorithm [Default: mmseqs-cluster]")
    parser_mmseqs.add_argument("-t", "--minimum_identity_threshold", type=float, default=50.0, help="Identity threshold [Default: 50.0]")
    parser_mmseqs.add_argument("--minimum_coverage_threshold", type=float, default=0.8, help="Coverage threshold [Default: 0.8]")
    parser_mmseqs.add_argument("--mmseqs2_options", type=str, default="", help="Additional MMseqs2 options")
    
    # Clustering arguments
    parser_clustering = parser.add_argument_group('Clustering arguments')
    parser_clustering.add_argument("--cluster_suffix", type=str, default="", help="Cluster suffix [Default: '']")
    parser_clustering.add_argument("--cluster_prefix_zfill", type=int, default=0, help="Prefix zfill [Default: 0]")
    parser_clustering.add_argument("--cluster_label_mode", type=str, default="md5", choices=["numeric", "random", "pseudo-random", "md5", "nodes"], help="Label mode [Default: md5]")
    parser_clustering.add_argument("--no_singletons", action="store_true", help="Exclude singletons")
    
    # Set the function to call
    parser.set_defaults(func=run)
    
    return parser

def run(args):
    """Execute cluster-proteins-from-pangenomes command"""
    
    # Validate
    if not args.pangenome_manifest and not (args.genome_manifest and args.genome_clusters):
        logger.error("Must provide either --pangenome_manifest OR both --genome_manifest and --genome_clusters")
        sys.exit(1)
    if args.pangenome_manifest and (args.genome_manifest or args.genome_clusters):
        logger.error("Cannot use --pangenome_manifest with --genome_manifest/--genome_clusters")
        sys.exit(1)
    
    if args.n_jobs == -1:
        args.n_jobs = cpu_count()
    
    # Setup directories with command-specific subdirectory
    directories = setup_directories(args.output_directory, subdirectory="pangenome_protein_clustering")
    
    # Create additional output subdirectories for comprehensive output
    os.makedirs(os.path.join(directories["output"], "pangenome_tables"), exist_ok=True)
    
    # Setup logger
    setup_logger(directories["log"], "cluster_proteins_from_pangenomes.log")
    
    # Print info
    logger.info("="*80)
    logger.info("pangenomium cluster-proteins-from-pangenomes")
    logger.info("="*80)
    print_header(version=__version__, n_jobs=args.n_jobs)
    
    # Parse input
    logger.info("-"*80)
    logger.info("Parsing input")
    logger.info("-"*80)
    pangenome_to_proteins = parse_input(args, directories)
    logger.info(f"Pangenomes: {len(pangenome_to_proteins)}")
    logger.info("")
    
    # Process each pangenome
    all_cluster_files = []
    algorithm = "easy-cluster" if args.algorithm == "mmseqs-cluster" else "easy-linclust"
    
    for pangenome_id, protein_data in tqdm(pangenome_to_proteins.items(), desc="Processing pangenomes", unit=" pangenomes"):
        protein_fasta = protein_data["fasta"]
        
        logger.info("="*80)
        logger.info(f"Clustering pangenome: {pangenome_id}")
        logger.info("="*80)
        
        # Run MMseqs2
        mmseqs_prefix = os.path.join(directories["intermediate"], f"{pangenome_id}_mmseqs2")
        edgelist_file = os.path.join(directories["intermediate"], f"{pangenome_id}_edgelist.tsv")
        
        cmd = [
            "mmseqs", algorithm,
            protein_fasta,
            mmseqs_prefix,
            directories["tmp"],
            "--threads", str(args.n_jobs),
            "--min-seq-id", str(args.minimum_identity_threshold / 100),
            "-c", str(args.minimum_coverage_threshold),
            "--cov-mode", "1",
        ]
        
        if args.mmseqs2_options:
            cmd.append(args.mmseqs2_options)
        
        cmd.extend([
            "&&",
            "mv", f"{mmseqs_prefix}_cluster.tsv", edgelist_file,
            "&&",
            "rm", "-rf",
            f"{mmseqs_prefix}_all_seqs.fasta",
            f"{mmseqs_prefix}_rep_seq.fasta",
        ])
        
        step = RunShellCommand(
            command=cmd,
            name=f"mmseqs2_{pangenome_id}",
            validate_output_filepaths=[edgelist_file]
        ).run()
        step.check_status()
        
        # Compile clusters with pangenome ID as prefix
        # This creates cluster IDs like: PanG-abc123__PanOG-xyz789
        cluster_file = os.path.join(directories["output"], f"{pangenome_id}_protein_clusters.tsv.gz")
        
        cmd = [
            "edgelist-to-clusters.py",
            "-i", edgelist_file,
            "-o", cluster_file,
            "--cluster_prefix", f"{pangenome_id}__PanOG-",  # Use pangenome ID + PanOG- as prefix
            "--cluster_prefix_zfill", str(args.cluster_prefix_zfill),
            "--cluster_label_mode", args.cluster_label_mode,
        ]
        
        # Only add cluster_suffix if non-empty
        if args.cluster_suffix:
            cmd.extend(["--cluster_suffix", args.cluster_suffix])
        
        if args.no_singletons:
            cmd.append("--no_singletons")
        
        step = RunShellCommand(
            command=cmd,
            name=f"compile_{pangenome_id}",
            validate_output_filepaths=[cluster_file]
        ).run()
        step.check_status()
        
        all_cluster_files.append(cluster_file)
        logger.info("")
    
    # Concatenate all cluster files
    if len(all_cluster_files) > 0:
        logger.info("="*80)
        logger.info("Concatenating cluster files")
        logger.info("="*80)
        final_output = os.path.join(directories["output"], "protein_clusters.tsv.gz")
        
        cmd = ["cat"] + all_cluster_files + [">", final_output]
        
        step = RunShellCommand(
            command=cmd,
            name="concatenate",
            validate_output_filepaths=[final_output]
        ).run()
        step.check_status()
        logger.info("")
    
    logger.info("="*80)
    logger.info("Complete")
    logger.info("="*80)
    return 0
